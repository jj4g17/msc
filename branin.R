branin <- function(XX, a=1, b=5.1/(4*pi^2), c=5/pi, r=6, s=10, t=1/(8*pi))
{
  x1 <- XX[1]; x2 <- XX[2]
  
  term1 <- a * (x2 - b*x1^2 + c*x1 - r)^2
  term2 <- s*(1-t)*cos(x1)
  
  y <- term1 + term2 + s
  return(y)
}

branin.analyse <- function(noise = 0, opt = "EI", 
                           design = list(lhs = optimumLHS, n0 = 6, rep = F), 
                           model = list(nugget = 0, MLE = F, update = F), 
                           sim = list(nruns=100, max=20, verb=1), f = branin) {
  
  performance <- data.frame(matrix(nrow=sim$nruns, ncol=7))
  colnames(performance) <- c("no_runs", "avg_opt_x1", "avg_opt_x2", "min_samp_Y", "min_mu", "nugget", "run_time")
  
  for (run in 1 : sim$nruns) {
    performance[run,] <- branin.analyse.seed(run, noise, opt, design, model, sim$max, sim$verb)
  }
  print(colMeans(performance))
  plot.dens(performance)
  return(performance)
}

branin.analyse.seed <- function(seed, noise, opt, design, model, max, verb, f = branin) {
  start.time <- proc.time()[3]; eps <- sqrt(.Machine$double.eps)
  
  set.seed(seed); if (verb > 0) cat("seed", seed, "\n")
  X1 <- seq(-5, 10, length=50); X2 <- seq(0, 15, length=50)
  XX <- expand.grid(X1, X2)
  
  X <- design$lhs(n0, 2)
  if (design$rep) X <- rbind(X, X)
  X[,1] <- -5 + 15*X[,1]
  X[,2] <- 15*X[,2]
  y <- apply(X, 1, f) + rnorm(dim(X)[1], 0, noise)
  d <- darg(NULL, X)
  if (model$MLE) g <- garg(list(mle=TRUE), y)
  else g <- list(start=model$nugget + eps)
  gpi <- newGP(X, y, d$start, g$start, dK=T)
  
  counter <- 0
  while (counter < max) {
    
    p <- predGP(gpi, XX, lite=T, nonug=T)
    
    if (model$update) g$start <- jmleGP(gpi, drange=c(sqrt(.Machine$double.eps), 5000), grange=c(sqrt(.Machine$double.eps), 5000))$g
    fmin <- min(y)
    #if (noise != 0 ) fmin <- min(predGP(gpi, X, lite = T)$mean)
    
    if (opt == "EI") {
      dist <- fmin - p$mean
      suppressWarnings(sigma <- sqrt(p$s2))
      distn <- dist/sigma
      surface <- dist*pnorm(distn) + sigma*dnorm(distn)
      newIndex <- which.max(surface)
    }
    
    if (opt == "ALC") surface <- alcGP(gpi, XX)
    
    newIndex <- which.max(surface)
    X.new <- XX[newIndex,]
    Y.new <- apply(X.new, 1, f) + rnorm(1, 0, noise)
    
    updateGP(gpi, X.new, Y.new)
    X <- rbind(X, t(matrix(X.new)))
    y <- c(y, Y.new)
    counter <- counter + 1
  }
  
  deleteGP(gpi)
  return(c(counter, XX[which.min(p$mean),1], XX[which.min(p$mean),2], min(y), min(p$mean), g$start, proc.time()[3] - start.time))
}

plot.dens <- function(obj, bins = 50, labs = F, obj2=NULL) {
  true.mins <- data.frame(x1=c(-pi, pi, 9.425), x2=c(12.275, 2.275, 2.475))
  X <- data.frame(x1=obj$avg_opt_x1, x2=obj$avg_opt_x2)
  if (is.null(obj2)) {
    ggplot(data=X, aes(x1, x2)) +
      stat_bin2d(bins = bins) +
      {if(labs)stat_bin2d(bins = bins, geom = 'text', colour="black", 
                          size=3.5, aes(label=..count..), 
                          position=position_nudge(x=0.75, y=0.05))} +
      geom_point(data=true.mins, aes(x1, x2), col = 'red', size = 5, pch=13)
  } else {
    X2 <- data.frame(x1=obj2$avg_opt_x1, x2=obj2$avg_opt_x2) 
    
    p1 <- ggplot(data=X, aes(x1, x2)) +
          stat_bin2d(bins = bins) +
          {if(labs)stat_bin2d(bins = bins, geom = 'text', colour="black", 
                              size=3.5, aes(label=..count..), 
                              position=position_nudge(x=0.75, y=0.05))} +
          geom_point(data=true.mins, aes(x1, x2), col = 'red', size = 5, pch=13)
    
    p2 <- ggplot(data=X2, aes(x1, x2)) +
          stat_bin2d(bins = bins) +
          {if(labs)stat_bin2d(bins = bins, geom = 'text', colour="black", 
                              size=3.5, aes(label=..count..), 
                              position=position_nudge(x=0.75, y=0.05))} +
          geom_point(data=true.mins, aes(x1, x2), col = 'red', size = 5, pch=13)
    
    grid.arrange(p1, p2, ncol=2)
    
  }
    

}



branin.EI <- branin.analyse()
plot(branin.EI$avg_opt_x1, branin.EI$avg_opt_x2, xlim=c(-5, 10), ylim=c(0, 15), xlab="x1", ylab='x2')
points(true.mins, pch=19, col='red', cex=1)

# no noise no nugget
EI.1 <- branin.analyse(noise = 0, opt = "EI", 
                       design=list(lhs=optimumLHS, n0 = 6, rep = F), 
                       sim = list(nruns = 1000, max = 20, verb = 2))
plot.dens(EI.1)
colMeans(EI.1)

# no noise no nugget - n0 = 10
EI.2 <- branin.analyse(noise = 0, opt = "EI", 
                       design=list(lhs=optimumLHS, n0 = 10, rep = F), 
                       sim = list(nruns = 1000, max = 40, verb = 2))
plot.dens(EI.2, bins = 30, labs=T)

#some noise no nugget
EI.3 <- branin.analyse(noise = 5, opt = "EI", 
                       design=list(lhs=optimumLHS, n0 = 6, rep = F), 
                       sim = list(nruns = 1000, max = 20, verb = 2))
plot.dens(EI.3, bins = 60, labs=F)

EI.3.2 <- branin.analyse(noise = 5, opt = "EI", 
                       design=list(lhs=optimumLHS, n0 = 10, rep = F), 
                       sim = list(nruns = 1000, max = 40, verb = 2))
plot.dens(EI.3.2, bins = 60, labs=F)
plot.dens(EI.3, obj2=EI.3.2, bins = 30, labs = F)


#some noise - nugget = 5
EI.4 <- branin.analyse(noise = 5, opt = "EI", 
                       design=list(lhs=optimumLHS, n0 = 6, rep = F), 
                       model=list(nugget=5, MLE=F, update=F),
                       sim = list(nruns = 1000, max = 20, verb = 2))

EI.4.2 <- branin.analyse(noise = 5, opt = "EI", 
                       design=list(lhs=optimumLHS, n0 = 10, rep = F), 
                       model=list(nugget=5, MLE=F, update=F),
                       sim = list(nruns = 1000, max = 40, verb = 2))

plot.dens(EI.4, obj2=EI.4.2, bins = 30, labs = F)

# some noise - nugget = MLE
EI.5 <- branin.analyse(noise = 5, opt = "EI", 
                       design=list(lhs=optimumLHS, n0 = 6, rep = F), 
                       model=list(nugget=5, MLE=T, update=F),
                       sim = list(nruns = 1000, max = 20, verb = 2))

EI.5.2 <- branin.analyse(noise = 5, opt = "EI", 
                       design=list(lhs=optimumLHS, n0 = 10, rep = F), 
                       model=list(nugget=5, MLE=T, update=F),
                       sim = list(nruns = 1000, max = 40, verb = 2))

plot.dens(EI.5, obj2=EI.5.2, bins = 30, labs = F)

# some noise - nugget = MLE - replicated design (same n0)
EI.6 <- branin.analyse(noise = 5, opt = "EI", 
                       design=list(lhs=optimumLHS, n0 = 6, rep = T), 
                       model=list(nugget=5, MLE=T, update=F),
                       sim = list(nruns = 1000, max = 20, verb = 2))

EI.6.2 <- branin.analyse(noise = 5, opt = "EI", 
                       design=list(lhs=optimumLHS, n0 = 10, rep = T), 
                       model=list(nugget=5, MLE=T, update=F),
                       sim = list(nruns = 1000, max = 40, verb = 2))

plot.dens(EI.6, obj2=EI.6.2, bins = 30, labs = F)


# ALC

# no noise, no nugget
ALC.1 <- branin.analyse(noise = 0, opt = "ALC", 
                        design=list(lhs=optimumLHS, n0 = 6, rep = F), 
                        sim = list(nruns = 100, max = 20, verb = 2))
plot.dens(ALC.1, bins = 20)
