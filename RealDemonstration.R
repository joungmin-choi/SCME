library(mvtnorm)
library(binaryLogic)
workpath = "./Work/MethyProject/Entropy"

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

calSmat = function(p.vec, r.vec, n)
{
  library(mvtnorm)
  library(binaryLogic)
  q.vec = rep(NA, 2 ^ n)
  state.lst = as.binary(0 : (2 ^ n - 1), n = n)
  state.mat = t(matrix(as.integer(unlist(state.lst)), nrow = n))
  lp = length(p.vec)
  lr = length(r.vec)
  S.mat = matrix(NA, lr, lp)
  
  for (i in 1 : lr)
  {
    r = r.vec[i]
    for (j in 1 : lp)
    {
      p = p.vec[j]
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
        for (k in 1 : 2 ^ n)
        {
          lower.vec = lower.mat[k, ]
          upper.vec = upper.mat[k, ]
          q.vec[k] = pmvnorm(lower = lower.vec, upper = upper.vec, sigma = rho.mat)
        }
        S.mat[i, j] = -sum(q.vec * log2(q.vec))
      }
      if (r == 1) S.mat[i, j] = -((1 - p) * log2(1 - p) + p * log2(p))
      if (p == 0 || p == 1) S.mat[i, j] = 0
    }
  }
  return(S.mat)
}

filtpat = function(data.df, depth)
{
  nsegraw = nrow(data.df)
  chrraw.vec = data.df$Chr
  locraw.vec = data.df$`4-CpG`
  pattraw.lst = vector("list", nsegraw)
  valid.vec = rep(F, nsegraw)
  data.lst = strsplit(data.df$MethylPattern, split = "[:;]")
  for (i in 1 : nsegraw)
  {
    datai = matrix(data.lst[[i]], nrow = 2)
    patt.vec = rep(datai[1, ], as.numeric(datai[2, ]))
    patt.mat = matrix(patt.vec, ncol = 1)
    patt.unlst = as.integer(unlist(strsplit(patt.mat, split = "")))
    patt = t(matrix(patt.unlst, nrow = 4))
    if (nrow(patt) >= depth)
    {
      valid.vec[i] = T
      pattraw.lst[[i]] = patt
    }
  }
  chr.vec = chrraw.vec[valid.vec]
  loc.vec = locraw.vec[valid.vec]
  patt.lst = pattraw.lst[valid.vec]
  return(list(chr.vec = chr.vec, loc.vec = loc.vec, patt.lst = patt.lst))
}

analypat = function(patt.lst)
{
  library(binaryLogic)
  nsite = 4
  state.lst = as.binary(0 : (2 ^ nsite - 1), n = nsite)
  state.mat = t(matrix(as.integer(unlist(state.lst)), nrow = nsite))
  
  nseg = length(patt.lst)
  nread.vec = rep(NA, nseg)
  MML.vec = rep(NA, nseg)
  OME.vec = rep(NA, nseg)
  pval.vec = rep(NA, nseg)
  r.vec = rep(NA, nseg)
  for (i in 1 : nseg)
  {
    patt.mat = patt.lst[[i]]
    nread = nrow(patt.mat)
    nread.vec[i] = nread
    MML.vec[i] = mean(patt.mat)
    cor.mat = cor(patt.mat)
    r.vec[i] = mean(cor.mat[upper.tri((cor.mat))])
    # distr.raw = table(apply(patt.mat, 1, paste, collapse = ""))
    # pmf = distr.raw / sum(distr.raw)
    # OME.vec[i] = -sum(pmf * log2(pmf))
    phat.vec = rep(MML.vec[i], nsite) # homogeneous
    q.vec = rep(0, 2 ^ nsite)
    o.vec = rep(0, 2 ^ nsite)
    for (j in 1 : 2 ^ nsite)
    {
      state.vec = state.mat[j, ]
      q.vec[j] = prod(phat.vec ^ state.vec * (1 - phat.vec) ^ (1 - state.vec))
      o.vec[j] = sum(apply(patt.mat, 1, identical, state.vec))
    }
    ix.vec = which(o.vec != 0)
    f.vec = o.vec[ix.vec] / nread
    OME.vec[i] = -sum(f.vec * log2(f.vec))
    
    e.vec = nread * q.vec
    ix.vec = which(e.vec != 0)
    stat = sum((o.vec[ix.vec] - e.vec[ix.vec]) ^ 2 / e.vec[ix.vec])
    # df = length(ix.vec) - 1 - sum(phat.vec != 1 & phat.vec != 0) # nonhomo
    df = length(ix.vec) - 1 - 1 # homo
    if (df > 0)
    {
      pval.vec[i] = pchisq(stat, df, lower.tail = F)
    }
  }
  return(list(nread.vec = nread.vec, MML.vec = MML.vec, OME.vec = OME.vec, pval.vec = pval.vec, r.vec = r.vec))
}

mergepat = function(lst1, lst2)
{
  loc1.df = data.frame(lst1$chr.vec, lst1$loc.vec)
  loc2.df = data.frame(lst2$chr.vec, lst2$loc.vec)
  nseg1 = nrow(loc1.df)
  nseg2 = nrow(loc2.df)
  tmp1.vec = rep(NA, nseg1)
  tmp2.vec = rep(NA, nseg2)
  for (i in 1 : nseg1)
  {
    tmp1.vec[i] = paste(loc1.df[i, ], collapse = "|")
  }
  for (i in 1 : nseg2)
  {
    tmp2.vec[i] = paste(loc2.df[i, ], collapse = "|")
  }
  ix1.vec = match(tmp1.vec, tmp2.vec)
  ix2.vec = match(tmp2.vec, tmp1.vec)
  if (sum(!is.na(ix1.vec)) != sum(!is.na(ix2.vec)))
  {
    stop("Something is wrong!")
  }
  nseg = sum(!is.na(ix1.vec))
  chr.vec = lst1$chr.vec[which(!is.na(ix1.vec))]
  loc.vec = lst1$loc.vec[which(!is.na(ix1.vec))]
  patt.lst = vector("list", nseg)
  pval.dmr = rep(NA, nseg)
  for (i in 1 : nseg)
  {
    ix1 = which(lst1$chr.vec == chr.vec[i] & lst1$loc.vec == loc.vec[i])
    ix2 = which(lst2$chr.vec == chr.vec[i] & lst2$loc.vec == loc.vec[i])
    mat1 = lst1$patt.lst[[ix1]]
    mat2 = lst2$patt.lst[[ix2]]
    patt.lst[[i]] = rbind(mat1, mat2)
    tab = rbind(c(sum(mat1), prod(dim(mat1)) - sum(mat1)), c(sum(mat2), prod(dim(mat2)) - sum(mat2)))
    pval.dmr[i] = chisq.test(tab)$p.value
  }
  return(list(chr.vec = chr.vec, loc.vec = loc.vec, patt.lst = patt.lst, pval.dmr.vec = pval.dmr))
}

simupat = function(nread.vec, MML.vec, r.vec)
{
  nseg = length(nread.vec)
  nsite = 4
  state.lst = as.binary(0 : (2 ^ nsite - 1), n = nsite)
  state.mat = t(matrix(as.integer(unlist(state.lst)), nrow = nsite))
  sMML.vec = rep(NA, nseg)
  OME.vec = rep(NA, nseg)
  pval.vec = rep(NA, nseg)
  set.seed(10)
  for (i in 1 : nseg)
  {
    nread = nread.vec[i]
    p = MML.vec[i]
    r = r.vec[i]
    if (!is.na(r) & r > 0 & !(all.equal(r, 1) == T))
    {
      r.mat = matrix(r, nsite, nsite)
      diag(r.mat) = 1
      p.vec = rep(p, nsite)
      X.mat = rCB(nread, p.vec, r.mat)
    } else if (!is.na(r) & (all.equal(r, 1) == T))
    {
      X = rbinom(nread, 1, p)
      X.mat = matrix(rep(X, nsite), ncol = nsite)
    } else
    {
      X.mat = matrix(NA, nread, nsite)
      for (j in 1 : nsite)
      {
        X.mat[, j] = rbinom(nread, 1, p)
      }
    }
    phat = mean(X.mat)
    sMML.vec[i] = phat
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
    pval.vec[i] = pchisq(stat, length(ix.vec) - 1 - 1, lower.tail = F)
    ix.vec = which(o.vec != 0)
    f.vec = o.vec[ix.vec] / nread
    OME.vec[i] = -sum(f.vec * log2(f.vec))
  }
  return(list(MML.vec = sMML.vec, OME.vec = OME.vec, pval.vec = pval.vec))
}

fn1 = paste(workpath, "/6week_female_neun+/SRR921832_trimmed_bismark.CpG.4CG.cnv.10X.txt", sep = "")
fn2 = paste(workpath, "/6week_female_neun-/SRR921839_trimmed_bismark.CpG.4CG.cnv.10X.txt", sep = "")
data1.tmp = read.table(fn1, header = T, sep = "\t", stringsAsFactors = F)
data1.df = data1.tmp[, c(1, 2, 9)]
colnames(data1.df) = c("Chr", "4-CpG", "MethylPattern")
data2.tmp = read.table(fn2, header = T, sep = "\t", stringsAsFactors = F)
data2.df = data2.tmp[, c(1, 2, 9)]
colnames(data2.df) = c("Chr", "4-CpG", "MethylPattern")

res1.1 = filtpat(data1.df, 20)
res1.2 = analypat(res1.1$patt.lst)
nread1.vec = res1.2$nread.vec
MML1.vec = res1.2$MML.vec
r1.vec = res1.2$r.vec
OME1.vec = res1.2$OME.vec
pval1.vec = res1.2$pval.vec
ix1.sig.vec = which(pval1.vec < 0.05)
ix1.nsig.vec = which(pval1.vec >= 0.05)
OME1.r0.7.vec = calSmat(MML1.vec, 0.7, 4)
ix1.bipo.vec = which(OME1.vec < OME1.r0.7.vec & MML1.vec < 0.8 & MML1.vec > 0.2)
print(paste("Total # of segments in NeuN+: ", length(pval1.vec), sep = "")) # 1012
print(paste("# of bipolar segments in NeuN+: ", length(ix1.bipo.vec), sep = "")) # 536
res2.1 = filtpat(data2.df, 20)
res2.2 = analypat(res2.1$patt.lst)
nread2.vec = res2.2$nread.vec
MML2.vec = res2.2$MML.vec
r2.vec = res2.2$r.vec
OME2.vec = res2.2$OME.vec
pval2.vec = res2.2$pval.vec
ix2.sig.vec = which(pval2.vec < 0.05)
ix2.nsig.vec = which(pval2.vec >= 0.05)
OME2.r0.7.vec = calSmat(MML2.vec, 0.7, 4)
ix2.bipo.vec = which(OME2.vec < OME2.r0.7.vec & MML2.vec < 0.8 & MML2.vec > 0.2)
print(paste("Total # of segments in NeuN-: ", length(pval2.vec), sep = "")) # 976
print(paste("# of bipolar segments in NeuN-: ", length(ix2.bipo.vec), sep = "")) # 527

res3.1 = mergepat(res1.1, res2.1)
ix3.dmr.vec = which(res3.1$pval.dmr.vec < 0.05) # 191
ix3.nondmr.vec = which(res3.1$pval.dmr.vec >= 0.05) # 775
sum(ix3.dmr.vec %in% ix3.bipo.vec) # 119

res3.2 = analypat(res3.1$patt.lst)
nread3.vec = res3.2$nread.vec
MML3.vec = res3.2$MML.vec
r3.vec = res3.2$r.vec
OME3.vec = res3.2$OME.vec
pval3.vec = res3.2$pval.vec
ix3.sig.vec = which(pval3.vec < 0.05)
ix3.nsig.vec = which(pval3.vec >= 0.05)
OME3.r0.7.vec = calSmat(MML3.vec, 0.7, 4)
ix3.bipo.vec = which(OME3.vec < OME3.r0.7.vec & MML3.vec < 0.8 & MML3.vec > 0.2)
ix3.nonbipo.vec = setdiff(1 : length(OME3.vec), ix3.bipo.vec)
tab = cbind(c(length(intersect(ix3.dmr.vec, ix3.bipo.vec)), length(intersect(ix3.dmr.vec, ix3.nonbipo.vec))), c(length(intersect(ix3.nondmr.vec, ix3.bipo.vec)), length(intersect(ix3.nondmr.vec, ix3.nonbipo.vec))))
fisher.test(tab)
