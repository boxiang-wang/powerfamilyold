hsvm <- function(v, varlist) {
  delta = varlist$delta
  r <- v[1]
  if (r > 1) 
    dl <- 0 else if (r <= (1 - delta)) 
      dl <- (1 - r - delta/2) else dl <- (r - 1)^2/delta/2
  dl
}

dhsvm <- function(v, varlist) {
  delta = varlist$delta
  r <- v[1]
  if (r > 1) 
    dl <- 0 else if (r <= (1 - delta)) 
      dl <- -1 else dl <- (r - 1) / delta
  dl
}

power <- function(v, varlist) {
  qv = varlist$qv
  decib = qv / (qv + 1)
  r <- v[1]
  dl = ifelse(r > decib, r ^ (-qv) * (qv ^ qv) / ((qv + 1) ^ (qv + 1)), 1 - r)
  dl
}

dpower <- function(v, varlist) {
  qv = varlist$qv
  decib = qv / (qv + 1)
  r <- v[1]
  dl = ifelse(r > decib, (-1) * r ^ (-qv - 1) * decib ^ (qv + 1), -1)
}


margin <- function(b0, beta, y, x, loss = c("hsvm", "power"), delta=2, qv=2) {
  loss <- match.arg(loss)
  nobs <- nrow(x)
  b0MAT <- matrix(rep(b0, nobs), nrow = nobs, byrow = TRUE)
  link <- x %*% beta + b0MAT
  if (loss %in% c("hsvm", "power")) {
    r <- y * link
  } else r <- y - link
  fun <- paste("d", loss, sep = "")
  varlist = list(delta=delta, qv=qv)
  dMat <- apply(r, c(1, 2), eval(fun), varlist) #dMat1
  if (loss %in% c("hsvm", "power")) {
    yxdMat <- t(x) %*% (dMat * y)/nobs
  } else yxdMat <- t(x) %*% dMat/nobs
  yxdMat
}
# l is the number of lambda
# p is nvars, n is nobs
# dim(dMat) = n by l
# dim(yxdMat) = p by l


#################################################################################
############ KKT, no printing, only give the violation percentage ################
#################################################################################

KKTnp = function(b0, beta, y, x, lambda, lambda2, thr, 
               loss = c("hsvm", "power"), delta=2, qv=2) {
  loss = match.arg(loss)
  dl = margin(b0, beta, y, x, loss=loss, delta=delta, qv=qv)
  ctr = 0
  for (l in 1:length(lambda)) {
    p = nrow(beta)
    for(j in 1:p)
    {
      if(beta[j,l]==0)
      {
        BB = abs(dl[j,l]) - lambda[l]
        if (BB > thr) 
        {
          ctr <- ctr + 1
        }
      } else{
        AA = dl[j,l] + lambda[l] * sign(beta[j,l]) + lambda2 * beta[j,l]
        if (abs(AA) >= thr)
        {
          ctr <- ctr + 1
        }
        
      }
    }
  }
  ccounts = length(lambda) * p
  return(ctr/ccounts*100)
}


#################################################################################
############ KKT, printing violations  ################
#################################################################################


KKTp = function(b0, beta, y, x, lambda, lambda2, thr, 
                loss = c("hsvm", "power"), delta=2, qv=2) {
  loss = match.arg(loss)
  dl = margin(b0, beta, y, x, loss=loss, delta=delta, qv=qv)
  count0 = 0
  ctr0 = 0
  ctrn0 = 0
  for (l in 1:length(lambda)) {
    p = nrow(beta)
    for(j in 1:p)
    {
      if(beta[j,l]==0)
      {
        BB = abs(dl[j,l]) - lambda[l]
        count0 = count0 + 1
        if (BB > thr) 
        {
          cat("violate at b = 0", BB, "\n")
          ctr0 <- ctr0 + 1
        }
      } else{
        AA = dl[j,l] + lambda[l] * sign(beta[j,l]) + lambda2 * beta[j,l]
        if (abs(AA) >= thr)
        {
          cat("violate at b != 0", abs(AA), "\n")
          ctrn0 <- ctrn0 + 1
        }
        
      }
    }
  }
  ctr = ctrn0 + ctr0
  ccounts = length(lambda) * p
  countn0 = ccounts - count0
  cat("# of checkings is ", ccounts, ".\n", sep="")
  cat("# of violations for zero beta is ", ctr0, ".\n", sep="")
  cat("% of violations for zero beta is ", ctr0/count0*100, "%", ".\n", sep="")
  cat("# of violations for non-zero beta is ", ctrn0, ".\n", sep="")
  cat("% of violations for non-zero beta is ", ctrn0/countn0*100, "%", ".\n", sep="")
  cat("# of total violations is ", ctr, ".\n", sep="")
  cat("% of total violations is ", ctr/ccounts*100, "%", ".\n", sep="")
  return(ctr/ccounts*100)
} 

#################################################################################
############ KKT, call KKTnp or KKTp by print.out ################
#################################################################################


KKT = function(b0, beta, y, x, lambda, lambda2, thr, 
               loss = c("hsvm", "power"), delta=2, qv=2, print.out = F)
{
  if(print.out == F)
  {
    KKTnp(b0=b0, beta=beta, y=y, x=x, lambda=lambda, 
                    lambda2=lambda2, thr=thr, 
                    loss = loss, delta=delta, qv=qv)
  } else
  {
    KKTp(b0=b0, beta=beta, y=y, x=x, lambda=lambda, 
                    lambda2=lambda2, thr=thr,  
                    loss = loss, delta=2, qv=2)
  }
}

#################################################################################
############ Fit model and check KKT conditions ################
#################################################################################


KKTperctg = function(dat, lambda2, qv, eps, thr)
{
  if(eps > 1) eps = 10 ^ (-eps)
  if(thr > 1) thr = 10 ^ (-thr)
  
  m.temp = GCDpower(x=dat$x, y=dat$y,
                       lambda2=lambda2, qv=qv, method="power",eps=eps, standardize=F)
  
  KKT(m.temp$b0, m.temp$beta, dat$y, dat$x, m.temp$lambda, lambda2=lambda2, thr=thr, 
      qv=qv, loss = c("power"), print.out=F)
}

#################################################################################
############ Summarize KKT condition checking tables ################
#################################################################################


KKTtb = function(dat, lambda2, qv, nm, eps.list=c(6:10), thr.list=c(2:5),
                 file.loc = "D:\\GitHub\\powerfamily\\Outputs\\KKTrda\\")
{
  perct.tb = matrix(NA, length(thr.list), length(eps.list))
  colnames(perct.tb) = paste("1e-", eps.list, sep="")
  rownames(perct.tb) = paste("1e-", thr.list, sep="")
  for(i in 1:length(eps.list))
  {
    for(j in 1:length(thr.list))
    {
      print(c(i, j))
      perct.tb[j, i] = KKTperctg(dat, lambda2=lambda2, qv=qv, eps=eps.list[i], 
                                 thr=thr.list[j])
    }
  }
  save("perct.tb", 
       file=paste(file.loc, nm, ".rda", sep="")) 
  return(perct.tb)
}
