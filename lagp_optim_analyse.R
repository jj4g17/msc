library(laGP); library(lhs); library(minimaxdesign); library(MaxPro); library(ggplot2); library(gridExtra);
eps <- sqrt(.Machine$double.eps); 

# design:
# LHS function | n0 | replication

# model:
# nugget | MLE | update |
# 

# sim:
# nruns | max | verb
lagp.optim.analyse <- function(f, xlim, noise=0, opt="EI",
                               design=list(lhs=optimumLHS, n0=6, rep=F),
                               model=list(nugget=0, MLE=F, update=F), sim=list(nruns=100, max=20, verb=1)) {
  
  performance <- data.frame(matrix(nrow=sim$nruns, ncol=6))
  colnames(performance) <- c("no_runs", "avg_opt_x", "dist_true_X", "min_samp_Y", "min_mu", "run_time")
  
  for (run in 1 : sim$nruns) performance[run,] <- lagp.optim.analyse.seed(run, f, xlim, noise, opt, design, model, sim$max, sim$verb)
  
  XX <- seq(xlim[1], xlim[2], length=1000)
  performance$`dist_true_X` <- abs(performance$`dist_true_X` - XX[which.min(f(XX))])
  
  return(performance)
}

lagp.optim.analyse.seed <- function(seed=1, f, xlim, noise=0, opt="EI", design=list(lhs=optimumLHS, n0=6, rep=F),
                                    model=list(nugget=0, MLE=F, update=F), max = 20, verb = 1) 
{
  start.time <- proc.time()[3]
  set.seed(seed); if (verb > 0) cat("Seed", seed, "\n")
  XX <- matrix(seq(xlim[1], xlim[2], length=1000))
  X <- matrix(design$lhs(design$n0, 1)); X <- xlim[1] + (xlim[2]-xlim[1])*X
  if (design$rep) X <- c(X, X)
  y <- f(X) + rnorm(length(X), 0, noise)
  d <- darg(NULL, X)
  g <- list(start=model$nugget + eps)
  if (model$MLE) g <- garg(list(mle=TRUE), y)
  gpi <- newGP(X, y, d$start, g$start, dK=T)
  
  counter <- 0
  while (counter < max) {
    if (verb > 2) print(counter)
    p <- predGP(gpi, XX, lite=T, nonug=T)
    if (model$update) jmleGP(gpi, drange=c(d$min, d$max), grange=c(eps, 1))
    fmin <- min(y)
    if (noise != 0 ) fmin <- min(predGP(gpi, matrix(X), lite = T)$mean)
    
    if (opt == "EI" | opt == "both") {
      dist <- fmin - p$mean
      suppressWarnings(sigma <- sqrt(p$s2))
      distn <- dist/sigma
      surface <- dist*pnorm(distn) + sigma*dnorm(distn)
      newIndex <- which.max(surface)
    }
    
    if (opt == "IECI" | opt == "both") surface <- 1-ieciGP(gpi, matrix(XX), fmin, nonug=T)
    
    if (opt == "ALC") surface <- alcGP(gpi, matrix(XX))
    
    newIndex <- which.max(surface)
    X.new <- XX[newIndex]
    Y.new <- f(X.new) + rnorm(1, 0, noise)
    
    if (verb > 2) ggplot.gp(gpi, X, y, f, X.new, Y.new, xlim, improvement=surface)
    
    updateGP(gpi, matrix(X.new), Y.new)
    X <- c(X, X.new)
    y <- c(y, Y.new)
    counter = counter + 1
  }
  
  if (verb > 1) ggplot.gp(gpi, X, y, f, X.new, Y.new, xlim, improvement=surface)
  
  deleteGP(gpi)
  
  return ( c(counter, XX[which.min(p$mean)], XX[which.min(p$mean)], min(y), min(p$mean), proc.time()[3] - start.time) )
}

ggplot.gp <- function(gpi, x, y, f, x.new, y.new, xlim, ylim, improvement) {
  XX <- seq(xlim[1], xlim[2], length=1000)
  yy <- f(XX)
  p <- predGP(gpi, matrix(XX), lite = T, nonug = T)
  df <- data.frame(XX=XX, mean=p$mean, s2 = p$s2, true=yy)
  df.pts <- data.frame(x=x, y=y)
  df.newpts <- data.frame(x.new=x.new, y.new=y.new)
  
  plot <- ggplot(df, aes(x=df$XX, y=df$mean)) +
    geom_line(color="#56B4E9") +
    geom_ribbon(aes(ymin=df$mean - 1.96*sqrt(df$s2), ymax = df$mean + 1.96*sqrt(df$s2)), alpha=0.1, linetype=2, color="#56B4E9") +
    geom_line(aes(x=df$XX, y=df$true), alpha=1, color="#CC79A7") +
    geom_point(data=df.pts, aes(x=df.pts$x, y=df.pts$y), size=2) +
    geom_point(data=df.newpts, aes(x=df.newpts$x.new, y=df.newpts$y.new), size=2, color='red') +
    labs(x='x', y='y')

  df.ei <- data.frame(x = XX, y=improvement)
  plot.EI <- ggplot(df.ei, aes(x=x, y=improvement)) +
    geom_line()
  grid.arrange(plot, plot.EI, ncol=1)
  readline()
  
}