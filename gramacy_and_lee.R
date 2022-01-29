grlee12 <- function(x)
{
  term1 <- sin(10*pi*x) / (2*x)
  term2 <- (x-1)^4
  
  y <- term1 + term2
  return(y)
}

plotDensity <- function(m, xlim, true_min) {
  ggplot(data=m, aes(x=avg_opt_x)) +
    geom_density(size=1) +
    geom_vline(xintercept=c(xlim[1], xlim[2], true_min), lty=2) +
    xlab('x') +
    xlim(xlim[1], xlim[2])
}

plotNuggetDensity <- function(m, xlim, true_min) {
  ggplot(data=m, aes(x=nugget)) +
    geom_density(size=1) +
    geom_vline(xintercept=c(xlim[1], xlim[2], true_min), lty=2) +
    xlab('g') +
    xlim(xlim[1], xlim[2])
} 

## true optimal - x = 0.549

## EI - no noise - no nugget
m.ei <- lagp.optim.analyse(f = grlee12, xlim = c(0.5, 2.5),
                              design = list(lhs = optimumLHS, n0=10, rep=F),
                              model = list(nugget = 0, MLE = F, update = F),
                              sim = list(nruns = 500, max = 20, verb = 1))

## EI - no noise - 0.5 nugget
m.ei.nug <- lagp.optim.analyse(f = grlee12, xlim = c(0.5, 2.5), noise = 0,
                           design = list(lhs = optimumLHS, n0=10, rep=F),
                           model = list(nugget = 0.5, MLE = F, update = F),
                           sim = list(nruns = 1000, max = 20, verb = 1))

ggplot(data=m.ei.nug, aes(x=avg_opt_x)) +
  geom_density(size=1) +
  geom_vline(xintercept=c(0.5,2.5, 0.548), lty=2) +
  xlab('x') +
  xlim(0.5, 2)

## EI - no noise - mle nugget
m.ei.nug <- lagp.optim.analyse(f = grlee12, xlim = c(0.5, 2.5), noise = 0,
                               design = list(lhs = optimumLHS, n0=10, rep=F),
                               model = list(nugget = 0.5, MLE = T, update = F),
                               sim = list(nruns = 1000, max = 20, verb = 1))

## EI - no noise - mle nugget - updates
m.ei.nug <- lagp.optim.analyse(f = grlee12, xlim = c(0.5, 2.5), noise = 0,
                               design = list(lhs = optimumLHS, n0=10, rep=F),
                               model = list(nugget = 0.5, MLE = T, update = T),
                               sim = list(nruns = 1000, max = 20, verb = 1))


## EI -  noise - nugget 0.5
m.ei.noise.nug <- lagp.optim.analyse(f = grlee12, xlim = c(0.5, 2.5), noise = 0.5,
                                     design = list(lhs = optimumLHS, n0=10, rep=F),
                                     model = list(nugget = 0.5, MLE = F, update = F),
                                     sim = list(nruns = 1000, max = 20, verb = 1))

m.ei.noise.nug.mle <- lagp.optim.analyse(f = grlee12, xlim = c(0.5, 2.5), noise = 0.5,
                                     design = list(lhs = optimumLHS, n0=10, rep=F),
                                     model = list(nugget = 0.5, MLE = F, update = F),
                                     sim = list(nruns = 100, max = 20, verb = 1))

m.ei.noise.nug.mle.rep <- lagp.optim.analyse(f = grlee12, xlim = c(0.5, 2.5), noise = 0.5,
                                     design = list(lhs = optimumLHS, n0=5, rep=T),
                                     model = list(nugget = 0.5, MLE = T, update = F),
                                     sim = list(nruns = 1000, max = 20, verb = 1))

m.ei.noise.nug.mle.rep.upd <- lagp.optim.analyse(f = grlee12, xlim = c(0.5, 2.5), noise = 0.5,
                                     design = list(lhs = optimumLHS, n0=5, rep=T),
                                     model = list(nugget = 0.5, MLE = T, update = T),
                                     sim = list(nruns = 1000, max = 20, verb = 1))

## EI - noise - nugget mle - 40 opts
m.ei.noise.nug.mle.more <- lagp.optim.analyse(f = grlee12, xlim = c(0.5, 2.5), noise = 0.5,
                                         design = list(lhs = optimumLHS, n0=10, rep=F),
                                         model = list(nugget = 0.5, MLE = T, update = F),
                                         sim = list(nruns = 1000, max = 40, verb = 1))


## IECI - no noise - no nugget mle
m.ieci <- lagp.optim.analyse(f = grlee12, xlim = c(0.5, 2.5), noise = 0, opt="IECI",
                               design = list(lhs = optimumLHS, n0=10, rep=F),
                               model = list(nugget = 0, MLE = F, update = F),
                               sim = list(nruns = 100, max = 20, verb = 1))

## IECI - noise - nugget 0.5
m.ieci.nug <- lagp.optim.analyse(f = grlee12, xlim = c(0.5, 2.5), noise = 0.5, opt="IECI",
                            design = list(lhs = optimumLHS, n0=10, rep=F),
                            model = list(nugget = 0.5, MLE = F, update = F),
                            sim = list(nruns = 100, max = 20, verb = 1))

## IECI - noise - nugget MLE
m.ieci.nug.mle <- lagp.optim.analyse(f = grlee12, xlim = c(0.5, 2.5), noise = 0.5, opt="IECI",
                                 design = list(lhs = optimumLHS, n0=10, rep=F),
                                 model = list(nugget = 0.5, MLE = T, update = F),
                                 sim = list(nruns = 100, max = 20, verb = 1))

## IECI - noise - nugget REPL
m.ieci.nug.mle.repl <- lagp.optim.analyse(f = grlee12, xlim = c(0.5, 2.5), noise = 0.5, opt="IECI",
                                     design = list(lhs = optimumLHS, n0=5, rep=T),
                                     model = list(nugget = 0.5, MLE = T, update = F),
                                     sim = list(nruns = 100, max = 20, verb = 1))

## IECI - noise - nugget REPL UPDATES
m.ieci.nug.mle.repl.upd <- lagp.optim.analyse(f = grlee12, xlim = c(0.5, 2.5), noise = 0.5, opt="IECI",
                                     design = list(lhs = optimumLHS, n0=5, rep=T),
                                     model = list(nugget = 0.5, MLE = T, update = T),
                                     sim = list(nruns = 100, max = 20, verb = 1))

m.ieci.nug.mle.upd <- lagp.optim.analyse(f = grlee12, xlim = c(0.5, 2.5), noise = 0.5, opt="IECI",
                                         design = list(lhs = optimumLHS, n0=10, rep=F),
                                         model = list(nugget = 0.5, MLE = T, update = T),
                                         sim = list(nruns = 100, max = 20, verb = 1))

## seed 29?
lagp.optim.analyse.seed(seed=29, xlim=c(0.5, 2.5), f=grlee12, opt="ALC", design = list(lhs = optimumLHS, n0=5, rep=T),
                        model = list(nugget = 0.5, MLE = T, update = T), verb=3)

## ALC - no noise - no nugget
m.alc <- lagp.optim.analyse(f = grlee12, xlim = c(0.5, 2.5), noise = 0, opt="ALC",
                             design = list(lhs = optimumLHS, n0=10, rep=F),
                             model = list(nugget = 0, MLE = F, update = F),
                             sim = list(nruns = 100, max = 20, verb = 2))
## ALC - noise - nugget 0.5
m.alc.nug <- lagp.optim.analyse(f = grlee12, xlim = c(0.5, 2.5), noise = 0.5, opt="ALC",
                            design = list(lhs = optimumLHS, n0=10, rep=F),
                            model = list(nugget = 0.5, MLE = F, update = F),
                            sim = list(nruns = 100, max = 20, verb = 2))

## ALC - noise - nugget mle
m.alc.nug.mle <- lagp.optim.analyse(f = grlee12, xlim = c(0.5, 2.5), noise = 0.5, opt="ALC",
                            design = list(lhs = optimumLHS, n0=10, rep=F),
                            model = list(nugget = 0.5, MLE = T, update = F),
                            sim = list(nruns = 100, max = 20, verb = 1))

## ALC - noise - nugget mle replicate
m.alc.nug.mle.repl <- lagp.optim.analyse(f = grlee12, xlim = c(0.5, 2.5), noise = 0.5, opt="ALC",
                            design = list(lhs = optimumLHS, n0=5, rep=T),
                            model = list(nugget = 0.5, MLE = T, update = F),
                            sim = list(nruns = 100, max = 20, verb = 2))

## ALC - noise - nugget mle - replicate -updates
m.alc.nug.mle.repl.upd <- lagp.optim.analyse(f = grlee12, xlim = c(0.5, 2.5), noise = 0.5, opt="ALC",
                                         design = list(lhs = optimumLHS, n0=5, rep=T),
                                         model = list(nugget = 0.5, MLE = T, update = T),
                                         sim = list(nruns = 100, max = 20, verb = 2))


