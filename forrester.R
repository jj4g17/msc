# full analysis of the forrester function
# http://www.sfu.ca/~ssurjano/forretal08.html
library(ggplot2)

forrester <- function(x) (6*x - 2)^2 * sin(12*x - 4)
XX <- seq(0, 1, length=10000)
plot(XX, forrester(XX), type='l', xlab="x", ylab='y')

# initial plot
ggplot(data=NULL, aes(x=XX, y=forrester(XX))) +
  geom_line(size=2) +
  xlab('x') +
  ylab('y')

# 6.2.1 : model 1 - EI, no nugget, Gauss kernel, deterministic design
evenDesign <- function(length, dim) seq(0, 1, length=length)
m1.lagp <- lagp.optim.analyse(f = forrester, xlim = c(0, 1),
                              design = list(lhs = optimumLHS, n0=6, rep=F),
                              model = list(nugget = 0, MLE = F, update = F),
                              sim = list(nruns = 1000, max = 20, verb = 1))
colMeans(m1.lagp)
var(m1.lagp$no_converge)
var(m1.lagp$avg_opt_x)

# EI, no nugget, Gauss kernel, noisy data sigma = 3
m1.ei.noisy <- lagp.optim.analyse(f = forrester, xlim=c(0, 1), noise = 3,
                                  design = list(lhs = evenDesign, n0=6, rep=F),
                                  model = list(nugget = 0, MLE = F, update = F),
                                  sim = list(nruns = 1000, max = 20, verb = 1))
colMeans(m1.ei.noisy)
var(m1.ei.noisy$no_converge)
var(m1.ei.noisy$avg_opt_x)

# EI, nugget = 3, sigma = 3
m1.ei.noisy.nugget <- lagp.optim.analyse(f = forrester, xlim=c(0, 1), noise = 3,
                                         design = list(lhs = optimumLHS, n0=6, rep=F),
                                         model = list(nugget = 3, MLE = F, update = F),
                                         sim = list(nruns = 1000, max = 20, verb = 1))

# EI, nugget = 3, sigma = 3, replicate
m1.ei.noisy.nugget.replicate <- lagp.optim.analyse(f = forrester, xlim=c(0, 1), noise = 3,
                                         design = list(lhs = optimumLHS, n0=3, rep=T),
                                         model = list(nugget = 3, MLE = F, update = F),
                                         sim = list(nruns = 1000, max = 20, verb = 2))

# EI, nugget = 3, sigma = 3, replicate, mle
m1.ei.noisy.nugget.replicate.mle <- lagp.optim.analyse(f = forrester, xlim=c(0, 1), noise = 3,
                                                   design = list(lhs = optimumLHS, n0=3, rep=T),
                                                   model = list(nugget = 3, MLE = T, update = F),
                                                   sim = list(nruns = 1000, max = 20, verb = 1))

# EI, nugget = 3, sigma = 3, replicate, mle update
m1.ei.noisy.nugget.replicate.mle.update <- lagp.optim.analyse(f = forrester, xlim=c(0, 1), noise = 3,
                                                       design = list(lhs = optimumLHS, n0=3, rep=T),
                                                       model = list(nugget = 3, MLE = T, update = T),
                                                       sim = list(nruns = 300, max = 20, verb = 2))
ggplot(data=m1.ei.noisy.nugget.replicate.mle.update, aes(x=avg_opt_x)) +
  geom_density(size=1) +
  geom_vline(xintercept=c(0,1, 0.757), lty=2) +
  xlab('x')


# EI, nugget = 3, sigma = 3, mle update
m1.ei.noisy.nugget.mle.update <- lagp.optim.analyse(f = forrester, xlim=c(0, 1), noise = 3,
                                                              design = list(lhs = optimumLHS, n0=6, rep=F),
                                                              model = list(nugget = 3, MLE = T, update = T),
                                                              sim = list(nruns = 300, max = 20, verb = 1))

ggplot(data=m1.ei.noisy.nugget.mle.update, aes(x=avg_opt_x)) +
  geom_density(size=1) +
  geom_vline(xintercept=c(0,1, 0.757), lty=2) +
  xlab('x')

# EI, nugget =3 , sigma = 3, mle update, spread design
m1.ei.noisy.nugget.mle.update <- lagp.optim.analyse(f = forrester, xlim=c(0, 1), noise = 3,
                                                    design = list(lhs = evenDesign, n0=6, rep=T),
                                                    model = list(nugget = 3, MLE = T, update = T),
                                                    sim = list(nruns = 50, max = 20, verb = 1))

# IECI, nugget = 0, sigma = 0, optimumLHS
m1.ieci <- lagp.optim.analyse(f = forrester, xlim=c(0, 1), noise = 0, opt="IECI",
                                                    design = list(lhs = optimumLHS, n0=6, rep=F),
                                                    model = list(nugget = 0, MLE = F, update = F),
                                                    sim = list(nruns = 100, max = 20, verb = 2))
ggplot(data=m1.ieci, aes(x=avg_opt_x)) +
  geom_density(size=1) +
  geom_vline(xintercept=c(0,1, 0.757), lty=2) +
  xlab('x')

# IECI, nugget = 3, sigma = 3
m1.ieci.nug3 <- lagp.optim.analyse(f = forrester, xlim=c(0, 1), noise = 3, opt="IECI",
                              design = list(lhs = optimumLHS, n0=6, rep=F),
                              model = list(nugget = 3, MLE = F, update = F),
                              sim = list(nruns = 100, max = 20, verb = 1))


# IECI nugget = MLE, sigma = 3
m1.ieci.nug3.mle <- lagp.optim.analyse(f = forrester, xlim=c(0, 1), noise = 3, opt="IECI",
                                   design = list(lhs = optimumLHS, n0=6, rep=F),
                                   model = list(nugget = 3, MLE = T, update = F),
                                   sim = list(nruns = 100, max = 40, verb = 1))

ggplot(data=m1.ieci.nug3.mle, aes(x=avg_opt_x)) +
  geom_density(size=1) +
  geom_vline(xintercept=c(0,1, 0.757), lty=2) +
  xlab('x')

# IECI nugget = mle, sigma = 3, double up
m1.ieci.nug3.mle.replicate <- lagp.optim.analyse(f = forrester, xlim=c(0, 1), noise = 3, opt="IECI",
                                       design = list(lhs = optimumLHS, n0=3, rep=T),
                                       model = list(nugget = 3, MLE = T, update = F),
                                       sim = list(nruns = 100, max = 20, verb = 2))
ggplot(data=m1.ieci.nug3.mle.replicate, aes(x=avg_opt_x)) +
  geom_density(size=1) +
  geom_vline(xintercept=c(0,1, 0.757), lty=2) +
  xlab('x')

# IECI nugget = mle, update MLE, double up
m1.ieci.nug3.mle.replicate.update <- lagp.optim.analyse(f = forrester, xlim=c(0, 1), noise = 3, opt="IECI",
                                                 design = list(lhs = optimumLHS, n0=6, rep=F),
                                                 model = list(nugget = 3, MLE = T, update = T),
                                                 sim = list(nruns = 100, max = 20, verb = 1))


# optimal model
m1.ieci.nug3.mle.replicate.update <- lagp.optim.analyse(f = forrester, xlim=c(0, 1), noise = 3, opt="IECI",
                                                        design = list(lhs = optimumLHS, n0=6, rep=T),
                                                        model = list(nugget = 3, MLE = T, update = T),
                                                        sim = list(nruns = 100, max = 20, verb = 1))

####################### ALC

# ALC - no noise - no nugget
alc <- lagp.optim.analyse(f = forrester, xlim=c(0, 1), opt="ALC",
                          sim=list(nruns=100, max=20, verb=1))

# ALC - noise - nugget = 3
alc.nug <- lagp.optim.analyse(f = forrester, xlim=c(0, 1), opt="ALC", noise = 3,
                          model = list(nugget = 3, MLE = F, update = F),
                          sim=list(nruns=100, max=20, verb=1))
colMeans(alc.nug)
ggplot(data=alc.nug, aes(x=avg_opt_x)) +
  geom_density(size=1) +
  geom_vline(xintercept=c(0,1, 0.757), lty=2) +
  xlab('x')
# ALC - noise - nugget = MLE
alc.nugMLE <- lagp.optim.analyse(f = forrester, xlim=c(0, 1), opt="ALC", noise = 3,
                          design = list(lhs = optimumLHS, n0= 6, rep = F),
                          model = list(nugget = 3, MLE = T, update = F),
                          sim=list(nruns=100, max=20, verb=1))
colMeans(alc.nugMLE)
ggplot(data=alc.nugMLE, aes(x=avg_opt_x)) +
  geom_density(size=1) +
  geom_vline(xintercept=c(0,1, 0.757), lty=2) +
  xlab('x')

# ALC - noise - MLE - double up
alc.nugMLEreplicate <- lagp.optim.analyse(f = forrester, xlim=c(0, 1), opt="ALC", noise = 3,
                          design = list(lhs = optimumLHS, n0= 3, rep = T),
                          model = list(nugget = 3, MLE = T, update = F),
                          sim=list(nruns=100, max=20, verb=1))
colMeans(alc.nugMLEreplicate)
ggplot(data=alc.nugMLEreplicate, aes(x=avg_opt_x)) +
  geom_density(size=1) +
  geom_vline(xintercept=c(0,1, 0.757), lty=2) +
  xlab('x')

# ALC - noise - MLE - update
alc.nugMLEUpdate <- lagp.optim.analyse(f = forrester, xlim=c(0, 1), opt="ALC", noise = 3,
                          design = list(lhs = optimumLHS, n0= 3, rep = T),
                          model = list(nugget = 3, MLE = T, update = T),
                          sim=list(nruns=100, max=20, verb=2))
ggplot(data=alc.nugMLEUpdate, aes(x=avg_opt_x)) +
  geom_density(size=1) +
  geom_vline(xintercept=c(0,1, 0.757), lty=2) +
  xlab('x')
colMeans(alc.nugMLEUpdate)
