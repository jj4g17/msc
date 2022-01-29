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
  
  performance <- data.frame(matrix(nrow=sim$nruns, ncol=8))
  colnames(performance) <- c("no_runs", "no_converge", "avg_opt_x", "dist_true_X", "min_samp_Y", "min_mu", "nugget", "run_time")
  
  for (run in 1 : sim$nruns) performance[run,] <- lagp.optim.analyse.seed(run, f, xlim, noise, opt, design, model, sim$max, sim$verb)
  
  XX <- seq(xlim[1], xlim[2], length=1000)
  performance$`dist_true_X` <- abs(performance$`dist_true_X` - XX[which.min(f(XX))])
  
  print(colMeans(performance))
  cat("\n Variance of final optimisers x: ", var(performance$avg_opt_x))
  
  return(performance)
}

lagp.optim.analyse.seed <- function(seed=1, f, xlim, noise=0, opt="EI", design=list(lhs=optimumLHS, n0=6, rep=F),
                                    model=list(nugget=0, MLE=F, update=F), max = 20, verb = 1) 
{
  start.time <- proc.time()[3]
  set.seed(seed); if (verb > 0) cat("Seed", seed, "\n")
  XX <- matrix(seq(xlim[1], xlim[2], length=1000))
  true.min.x <- XX[which.min(f(XX))]
  X <- matrix(design$lhs(design$n0, 1)); X <- xlim[1] + (xlim[2]-xlim[1])*X
  if (design$rep) X <- matrix(c(X, X))
  y <- f(X) + rnorm(length(X), 0, noise)
  d <- darg(NULL, X)
  
  if (model$MLE) g <- garg(list(mle=TRUE), y)
  else g <- list(start=model$nugget + eps)
  gpi <- newGP(X, y, d$start, g$start, dK=T)
  
  counter <- 0
  converge <- F
  convergePoint <- max
  while (counter < max) {
    if (verb > 2) print(counter)
    p <- predGP(gpi, XX, lite=T, nonug=T)
    if (model$update) {
      g$start <- jmleGP(gpi, grange=c(sqrt(.Machine$double.eps), 20))$g
    }
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
    
    if (verb > 2) ggplot.gp(gpi, X, y, f, X.new, Y.new, xlim, improvement=surface, plotSurface=T)
    
    updateGP(gpi, matrix(X.new), Y.new)
    X <- c(X, X.new)
    y <- c(y, Y.new)
    counter = counter + 1
    
    if (abs(XX[which.min(p$mean)] - true.min.x) < 0.001 & converge == F) {
      converge <- T
      convergePoint <- counter
    }
  }
  
  if (verb > 1) {
    ggplot.gp(gpi, X, y, f, X.new, Y.new, xlim, improvement=surface, plotSurface=F)
  }
  deleteGP(gpi)
  
  return ( c(counter, convergePoint, XX[which.min(p$mean)], XX[which.min(p$mean)], min(y), min(p$mean), g$start, proc.time()[3] - start.time) )
}

ggplot.gp <- function(gpi, x, y, f, x.new, y.new, xlim, ylim, improvement, plotSurface = F) {
  XX <- seq(xlim[1], xlim[2], length=1000)
  yy <- f(XX)
  p <- predGP(gpi, matrix(XX), lite = T, nonug = T)
  df <- data.frame(XX=XX, mean=p$mean, s2 = p$s2, true=yy)
  df.pts <- data.frame(x=x, y=y)
  df.newpts <- data.frame(x.new=x.new, y.new=y.new)
  
  plot.gp <- ggplot(df, aes(x=XX, y=mean)) +
    geom_line(color="#56B4E9") +
    geom_ribbon(data=df, aes(ymin=mean - 1.96*sqrt(s2), ymax = mean + 1.96*sqrt(s2)), alpha=0.1, linetype=2, color="#56B4E9") +
    geom_line(data=df, aes(x=XX, y=true), alpha=1, color="#CC79A7") +
    geom_point(data=df.pts, aes(x=x, y=y), size=2) +
    geom_point(data=df.newpts, aes(x=x.new, y=y.new), size=2, color='red') +
    labs(x='x', y='y')
  
  df.ei <- data.frame(x = XX, y=improvement)
  plot.EI <- ggplot(df.ei, aes(x=x, y=improvement)) +
    geom_line()
  
  if (plotSurface == T) {
    grid.arrange(plot.gp, plot.EI, ncol=1)
    readline()
  } 
  if (plotSurface == F) {
    grid.arrange(plot.gp, ncol=1, nrow=1)
    readline()
  }
}