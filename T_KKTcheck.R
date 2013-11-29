rm(list=ls(all=TRUE))
setwd("D:\\GitHub\\powerfamily")

dhsvm <- function(v, delta) {
  r <- v[1]
  if (r > 1) 
    dl <- 0 else if (r <= (1 - delta)) 
      dl <- (1 - r - delta/2) else dl <- (r - 1)^2/delta
  dl
}



margin <- function(b0, beta, y, x, delta, loss = c("ls", "logit", 
                                                   "sqsvm", "hsvm")) {
  loss <- match.arg(loss)
  nobs <- nrow(x)
  b0MAT <- matrix(rep(b0, nobs), nrow = nobs, byrow = TRUE)
  link <- x %*% beta + b0MAT
  if (loss %in% c("logit", "sqsvm", "hsvm")) {
    r <- y * link
  } else r <- y - link
  fun <- paste("d", loss, sep = "")
  dMat <- apply(r, c(1, 2), eval(fun), delta = delta) #dMat1
  if (loss %in% c("logit", "sqsvm", "hsvm")) {
    yxdMat <- t(x) %*% (dMat * y)/nobs
  } else yxdMat <- t(x) %*% dMat/nobs
  yxdMat
}
# l is the number of lambda
# p is nvars, n is nobs
# dim(dMat) = n by l
# dim(yxdMat) = p by l

KKT1 = function(b0, beta, y, x, lambda, thr, delta, loss = c("ls", 
                                                            "logit", "sqsvm", "hsvm")) {
  loss = match.arg(loss)
  B = matrix(NA, ncol = length(lambda))
  dl = margin(b0, beta, y, x, delta, loss)
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
          cat("violate at b = 0", BB, "\n")
          ctr <- ctr + 1
        }
      } else{
        AA = dl[j,l] + lambda[l] * sign(beta[j,l])
        if (abs(sum(AA)) >= thr)
        {
          cat("violate at b != 0", abs(sum(AA)), "\n")
          ctr <- ctr + 1
        }
        
      }
    }
  }
  cat("# of violations", ctr/length(lambda), "\n")
  return(ctr/length(lambda))
} 



library(gcdnet)
load("FHT.rda")
y = FHT$y
x = FHT$x



delta = as.double(2)
thr = 1e-04
b0=m$b0
beta=m$beta
lambda=m$lambda

source("GCDpower.R")
m <- gcdnet(x=FHT$x,y=FHT$y,
            #lambda=c(0.5,0.1),
            lambda2=0, delta=2,method="hhsvm")
plot(m, color=T)


margin(m$b0, m$beta, FHT$y, FHT$x, delta, loss = c("hsvm") )
KKT1(m$b0, m$beta, FHT$y, FHT$x, m$lambda, thr, delta, loss = c("hsvm"))
