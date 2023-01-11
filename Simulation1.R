# Simulation 1: evaluate the chi-square test for spatial correlation
library(binaryLogic)

# Type-I error
nsimu = 1e4
pval.vec = rep(NA, nsimu)
nsite.vec = c(2, 3, 4)
lns = length(nsite.vec)
nread.vec = c(20, 40, 60, 80, 100)
lnr = length(nread.vec)
typeI.mat = matrix(NA, lns, lnr)
set.seed(10)
system.time({
for (a in 1 : lns)
{
  nsite = nsite.vec[a]
  for (b in 1 : lnr)
  {
    nread = nread.vec[b]
    for (k in 1 : nsimu)
    {
      p.vec = runif(nsite, 0.3, 0.7)
      # p.vec = runif(nsite, 0.1, 0.9)
      X.mat = matrix(NA, nread, nsite)
      for (j in 1 : nsite)
      {
        X.mat[, j] = rbinom(nread, 1, p.vec[j])
      }
      
      phat.vec = colMeans(X.mat)
      state.lst = as.binary(0 : (2 ^ nsite - 1), n = nsite)
      state.mat = t(matrix(as.integer(unlist(state.lst)), nrow = nsite))
      q.vec = rep(NA, 2 ^ nsite)
      o.vec = rep(0, 2 ^ nsite)
      for (i in 1 : 2 ^ nsite)
      {
        state.vec = state.mat[i, ]
        q.vec[i] = prod(phat.vec ^ state.vec * (1 - phat.vec) ^ (1 - state.vec))
        o.vec[i] = sum(apply(X.mat, 1, identical, state.vec))
      }
      
      e.vec = nread * q.vec
      ix.vec = which(e.vec != 0)
      stat = sum((o.vec[ix.vec] - e.vec[ix.vec]) ^ 2 / e.vec[ix.vec])
      pval.vec[k] = pchisq(stat, length(ix.vec) - 1 - sum(phat.vec != 1 & phat.vec != 0), lower.tail = F)
    }
    typeI.mat[a, b] = mean(pval.vec < 0.05)
  }
}
})
typeI.mat # [0.0457, 0.0543]
# p \in (0, 1) run time: 844 secs
# [,1]   [,2]   [,3]   [,4]   [,5]
# [1,] 0.0795 0.0672 0.0648 0.0575 0.0633
# [2,] 0.0588 0.0529 0.0543 0.0559 0.0567
# [3,] 0.0693 0.0713 0.0711 0.0721 0.0696
# p \in (0.1, 0.9) run time: 444 secs
# [,1]   [,2]   [,3]   [,4]   [,5]
# [1,] 0.0486 0.0506 0.0474 0.0500 0.0516
# [2,] 0.0546 0.0445 0.0499 0.0496 0.0489
# [3,] 0.0655 0.0628 0.0595 0.0568 0.0493
# p \in (0.3, 0.7) run time: 391 secs
# [,1]   [,2]   [,3]   [,4]   [,5]
# [1,] 0.0489 0.0565 0.0505 0.0508 0.0503
# [2,] 0.0459 0.0467 0.0497 0.0475 0.0504
# [3,] 0.0455 0.0472 0.0502 0.0472 0.0520

# Power
library(mvtnorm)
library(binaryLogic)

posdef <- function(rho.mat)
{
  rho.mat.eig <- eigen(rho.mat)
  rho.mat.eva <- rho.mat.eig$values
  rho.mat.eve <- rho.mat.eig$vectors
  rho.mat.eva[rho.mat.eva < 0] <- 0
  rho.mat.rev <- rho.mat.eve %*% diag(rho.mat.eva) %*% t(rho.mat.eve)
  
  return(rho.mat.rev)
}

rCB <- function(n, p.vec, r.mat)
{
  # this function generates random samples from a multivariate cor-Bern random vector
  # Y = (Y1,...,Ym) where p(Yi = 1) = pi, and cor(Yi, Yj) = rij
  # note, the function uses a thresholding method (Emrich & Piedmonte)
  # note, use this function to simulate binary trait
  # input: 
  #         n: a scalar, the number of cor-Bern samples to generate, each sample is a length-m vector
  #         p.vec: a length-m vector, the Bernoulli mean
  #         r.mat: a m by m matrix, the correlation matrix
  # output:
  #         a n by m matrix of cor-Bern random samples
  
  m <- length(p.vec)
  if (nrow(r.mat) != m | ncol(r.mat) != m)
  {
    stop("wrong dimension!")
  }else
  {
    d.vec <- qnorm(p.vec, lower.tail = F)
    
    f <- function(rho, r, p1, p2)
    {
      d1 <- qnorm(p1, lower.tail = F)
      d2 <- qnorm(p2, lower.tail = F)
      s <- cbind(c(1, rho), c(rho, 1))
      E12 <- r * sqrt(p1 * (1 - p1) * p2 * (1 - p2)) + p1 * p2
      return(pmvnorm(lower = c(d1, d2), sigma = s) - E12)
    }
    
    rho.mat <- matrix(NA, m, m)
    diag(rho.mat) <- 1
    for (i in 1 : (m - 1))
    {
      for (j in (i + 1) : m)
      {
        r <- r.mat[i, j]
        rho <- uniroot(f, r = r, p1 = p.vec[i], p2 = p.vec[j], lower = -1, upper = 1, tol = 1e-9)$root
        rho.mat[i, j] <- rho
        rho.mat[j, i] <- rho
      }
    }
    rho.mat <- posdef(rho.mat)
    
    x <- rmvnorm(n, sigma = rho.mat)
    y <- matrix(as.integer(as.vector(x) >= rep(d.vec, each = n)), ncol = m)
    
    return(y)
  }
}

ar1 = function(n, rho)
{
  exponent = abs(matrix(1 : n - 1, nrow = n, ncol = n, byrow = T) - (1 : n - 1))
  return(rho ^ exponent)
}

r = 0.3 # 0.4
nsimu = 1e3
pval.vec = rep(NA, nsimu)
nsite.vec = c(2, 3, 4)
lns = length(nsite.vec)
nread.vec = c(20, 40, 60, 80, 100)
lnr = length(nread.vec)
power.mat = matrix(NA, lns, lnr)
set.seed(10)
system.time({
for (a in 1 : lns)
{
  nsite = nsite.vec[a]
  r.mat = ar1(nsite, r)
  for (b in 1 : lnr)
  {
    nread = nread.vec[b]
    for (k in 1 : nsimu)
    {
      p.vec = runif(nsite, 0.3, 0.7)
      X.mat = rCB(nread, p.vec, r.mat)

      phat.vec = colMeans(X.mat)
      state.lst = as.binary(0 : (2 ^ nsite - 1), n = nsite)
      state.mat = t(matrix(as.integer(unlist(state.lst)), nrow = nsite))
      q.vec = rep(NA, 2 ^ nsite)
      o.vec = rep(0, 2 ^ nsite)
      for (i in 1 : 2 ^ nsite)
      {
        state.vec = state.mat[i, ]
        q.vec[i] = prod(phat.vec ^ state.vec * (1 - phat.vec) ^ (1 - state.vec))
        o.vec[i] = sum(apply(X.mat, 1, identical, state.vec))
      }
      
      e.vec = nread * q.vec
      ix.vec = which(e.vec != 0)
      stat = sum((o.vec[ix.vec] - e.vec[ix.vec]) ^ 2 / e.vec[ix.vec])
      pval.vec[k] = pchisq(stat, length(ix.vec) - 1 - sum(phat.vec != 1 & phat.vec != 0), lower.tail = F)
    }
    power.mat[a, b] = mean(pval.vec < 0.05)
  }
}
})
power.mat
