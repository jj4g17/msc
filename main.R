
# LHS method comparison 6.1
library(lhs)

T <- 100; n0 <- 5
MM <- matrix(NA, nrow=3, ncol=4)
for (dim in 1:3) 
{
  M <- matrix(NA, nrow=T, ncol=4)
  for (t in 1:T) 
  {
    print(t)
    
    x.genetic <- geneticLHS(n0^dim, dim)
    x.maximin <- maximinLHS(n0^dim, dim)
    x.optimum <- optimumLHS(n0^dim, dim)
    x.random <-  randomLHS(n0^dim, dim)
    
    x.genetic.md <- min(dist(x.genetic))
    x.maximin.md <- min(dist(x.maximin))
    x.optimum.md <- min(dist(x.optimum))
    x.random.md <- min(dist(x.random))
    
    M[t,] <- c(x.random.md, x.maximin.md, x.optimum.md, x.genetic.md)
  }
  MM[dim,] <- colMeans(M)
}


