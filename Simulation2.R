# Simulation 2: validate the relation between OME and ML
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
      # Let E12 = E[Y1 Y2]. Since Cov(Y1, Y2) = E[Y1 Y2] - E[Y1]E[Y2]
      #                                       = r sqrt(var(Y1) var(Y2)),
      # E12 = r sqrt(var(Y1) var(Y2)) + E[Y1]E[Y2]
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

calSn = function(p, r, n)
{
  q.vec = rep(NA, 2 ^ n)
  state.lst = as.binary(0 : (2 ^ n - 1), n = n)
  state.mat = t(matrix(as.integer(unlist(state.lst)), nrow = n))
  if (r != 1 && p != 0 && p != 1)
  {
    delta = qnorm(p, lower.tail = F)
    f = function(rho, r, p)
    {
      d = qnorm(p, lower.tail = F)
      s <- cbind(c(1, rho), c(rho, 1))
      E12 <- r * p * (1 - p) + p ^ 2
      return(pmvnorm(lower = c(d, d), sigma = s) - E12)
    }
    rho = uniroot(f, r = r, p = p, lower = -1, upper = 1, tol = 1e-9)$root
    rho.mat = matrix(rho, n, n)
    diag(rho.mat) = 1
    lower.bd = c(-Inf, delta)
    lower.mat = matrix(lower.bd[state.mat + 1], ncol = n)
    upper.bd = c(delta, Inf)
    upper.mat = matrix(upper.bd[state.mat + 1], ncol = n)
    for (i in 1 : 2 ^ n)
    {
      lower.vec = lower.mat[i, ]
      upper.vec = upper.mat[i, ]
      q.vec[i] = pmvnorm(lower = lower.vec, upper = upper.vec, sigma = rho.mat)
    }
    S = -sum(q.vec * log2(q.vec))
  }
  if (r == 1) S = -((1 - p) * log2(1 - p) + p * log2(p))
  if (p == 0 || p == 1) S = 0
  return(S)
}

pp.vec = seq(0, 1, by = 0.05)
lp = length(pp.vec)
rr.vec = seq(0, 1, by = 0.1)
lr = length(rr.vec)
S.mat = matrix(NA, lr, lp)
for (i in 1 : lr)
{
  rr = rr.vec[i]
  for (j in 1 : lp)
  {
    pp = pp.vec[j]
    S.mat[i, j] = calSn(pp, rr, 4)
  }
}

nseg = 500
nullprop = 0.1
ixnull.vec = sample(nseg, nseg * nullprop)
ixalt.vec = setdiff(1 : nseg, ixnull.vec)
nsite = 4
nread.vec = c(20, 40, 60, 80, 100)
lnr = length(nread.vec)
state.lst = as.binary(0 : (2 ^ nsite - 1), n = nsite)
state.mat = t(matrix(as.integer(unlist(state.lst)), nrow = nsite))
nsig.mat = matrix(NA, lnr, 2)
title.vec = c("A", "B", "C", "D", "E")
x11()
par(mfrow = c(3, 2), cex.main = 2, cex.axis = 1.5, cex.lab = 1.5)
for (k in 1 : lnr)
{
  nread = nread.vec[k]
  MML.vec = rep(NA, nseg)
  OME.vec = rep(NA, nseg)
  pval.vec = rep(NA, nseg)
  set.seed(10)
  for (i in 1 : nseg)
  {
    p = runif(1, 0, 1)
    p.vec = rep(p, nsite)
    if (i %in% ixnull.vec)
    {
      X.mat = matrix(NA, nread, nsite)
      for (j in 1 : nsite)
      {
        X.mat[, j] = rbinom(nread, 1, p.vec[j])
      }
    } else
    {
      r = runif(1, 0, 1)
      r.mat = matrix(r, nsite, nsite)
      diag(r.mat) = 1
      X.mat = rCB(nread, p.vec, r.mat)
    }
    phat = mean(X.mat)
    MML.vec[i] = phat
    q.vec = rep(0, 2 ^ nsite)
    o.vec = rep(0, 2 ^ nsite)
    for (j in 1 : 2 ^ nsite)
    {
      state.vec = state.mat[j, ]
      q.vec[j] = prod(phat ^ state.vec * (1 - phat) ^ (1 - state.vec))
      o.vec[j] = sum(apply(X.mat, 1, identical, state.vec))
    }
    e.vec = nread * q.vec
    ix.vec = which(e.vec != 0)
    stat = sum((o.vec[ix.vec] - e.vec[ix.vec]) ^ 2 / e.vec[ix.vec])
    if (length(ix.vec) > 2)
    {
      pval.vec[i] = pchisq(stat, length(ix.vec) - 1 - 1, lower.tail = F)
    }
    ix.vec = which(o.vec != 0)
    f.vec = o.vec[ix.vec] / nread
    OME.vec[i] = -sum(f.vec * log2(f.vec))
  }
  nosigix.vec = which(pval.vec >= 0.05)
  plot(MML.vec[ixnull.vec], OME.vec[ixnull.vec], type = "p", pch = 20, col = "red", xlim = c(0, 1), ylim = c(0, 4), xlab = "MML", ylab = "OME", main = paste(title.vec[k], ": m=", nread, sep = ""))
  points(MML.vec[nosigix.vec], OME.vec[nosigix.vec], pch = 1, col = "blue")
  sigix.vec = which(pval.vec < 0.05)
  nsig.mat[k, 1] = length(sigix.vec)
  nsig.mat[k, 2] = length(nosigix.vec)
  points(MML.vec[sigix.vec], OME.vec[sigix.vec], pch = 20, col = "black")

  for (i in 2 : lr)
  {
    lines(pp.vec, S.mat[i, ], col = "black")
  }
  lines(pp.vec, S.mat[1, ], lwd = 2, col = "blue")
  text(0.75, 0.5, "r=1", cex = 1.3) # 1.5
  text(0.75, 3.75, "r=0", cex = 1.3) # 1.5
}
