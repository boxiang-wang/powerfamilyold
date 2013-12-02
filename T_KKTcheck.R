rm(list=ls(all=TRUE))
setwd("D:\\GitHub\\powerfamily")


#################################################################################
        ############ checking KKT conditions for HHSVM ################
#################################################################################
hsvm <- function(v, delta) {
  r <- v[1]
  if (r > 1) 
    dl <- 0 else if (r <= (1 - delta)) 
      dl <- (1 - r - delta/2) else dl <- (r - 1)^2/delta/2
  dl
}

dhsvm <- function(v, delta) {
  r <- v[1]
  if (r > 1) 
    dl <- 0 else if (r <= (1 - delta)) 
      dl <- -1 else dl <- (r - 1) / delta
  dl
}


margin <- function(b0, beta, y, x, delta, loss = c("hsvm")) {
  loss <- match.arg(loss)
  nobs <- nrow(x)
  b0MAT <- matrix(rep(b0, nobs), nrow = nobs, byrow = TRUE)
  link <- x %*% beta + b0MAT
  if (loss %in% c("hsvm")) {
    r <- y * link
  } else r <- y - link
  fun <- paste("d", loss, sep = "")
  dMat <- apply(r, c(1, 2), eval(fun), delta = delta) #dMat1
  if (loss %in% c("hsvm")) {
    yxdMat <- t(x) %*% (dMat * y)/nobs
  } else yxdMat <- t(x) %*% dMat/nobs
  yxdMat
}
# l is the number of lambda
# p is nvars, n is nobs
# dim(dMat) = n by l
# dim(yxdMat) = p by l

KKT1 = function(b0, beta, y, x, lambda, thr, delta, loss = c("hsvm")) {
  loss = match.arg(loss)
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
load("D_FHT.rda")
y = FHT$y
x = FHT$x

#x = matrix(c(1,0,-1,-1,1,0),3,2)
#y = c(-1,-1,1)


delta = as.double(2)
thr = 1e-04


#source("GCDpower.R")
m <- gcdnet(x=FHT$x,y=FHT$y,
            #lambda=c(0.1,0.01),
            lambda2=0, delta=2,method="hhsvm",eps=1e-10, standardize=F)
plot(m, color=T)

b0=m$b0
beta=m$beta
lambda=m$lambda


margin(m$b0, m$beta, FHT$y, FHT$x, delta=2, loss = c("hsvm") )
KKT1(m$b0, m$beta, FHT$y, FHT$x, m$lambda, thr=1e-03, delta=2, loss = c("hsvm"))


#################################################################################
############ checking KKT conditions for power family ################
#################################################################################
hsvm <- function(v, delta) {
  r <- v[1]
  if (r > 1) 
    dl <- 0 else if (r <= (1 - delta)) 
      dl <- (1 - r - delta/2) else dl <- (r - 1)^2/delta/2
  dl
}

dhsvm <- function(v, delta) {
  r <- v[1]
  if (r > 1) 
    dl <- 0 else if (r <= (1 - delta)) 
      dl <- -1 else dl <- (r - 1) / delta
  dl
}

power <- function(v, qv) {
  decib = qv / (qv + 1)
  r <- v[1]
  dl = ifelse(r > decib, r ^ (-qv) * (qv ^ qv) / ((qv + 1) ^ (qv + 1)), 1 - r)
  dl
}

dpower <- function(v, qv) {
  decib = qv / (qv + 1)
  r <- v[1]
  dl = ifelse(r > decib, (-1) * r ^ (-qv - 1) * decib ^ (qv + 1), -1)
}


margin <- function(b0, beta, y, x, delta, loss = c("hsvm", "power")) {
  loss <- match.arg(loss)
  nobs <- nrow(x)
  b0MAT <- matrix(rep(b0, nobs), nrow = nobs, byrow = TRUE)
  link <- x %*% beta + b0MAT
  if (loss %in% c("hsvm", "power")) {
    r <- y * link
  } else r <- y - link
  fun <- paste("d", loss, sep = "")
  dMat <- apply(r, c(1, 2), eval(fun), delta = delta) #dMat1
  if (loss %in% c("hsvm", "power")) {
    yxdMat <- t(x) %*% (dMat * y)/nobs
  } else yxdMat <- t(x) %*% dMat/nobs
  yxdMat
}
# l is the number of lambda
# p is nvars, n is nobs
# dim(dMat) = n by l
# dim(yxdMat) = p by l

KKT1 = function(b0, beta, y, x, lambda, lambda2, thr, delta, loss = c("hsvm")) {
  loss = match.arg(loss)
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
        AA = dl[j,l] + lambda[l] * sign(beta[j,l]) * lambda2 * beta[j,l]
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
load("D_FHT.rda")
y = FHT$y
x = FHT$x

#x = matrix(c(1,0,-1,-1,1,0),3,2)
#y = c(-1,-1,1)

#source("GCDpower.R")
m <- gcdnet(x=FHT$x,y=FHT$y,
            #lambda=c(0.1,0.01),
            lambda2=1, delta=2,method="hhsvm",eps=1e-10, standardize=F)
plot(m, color=T)

b0=m$b0
beta=m$beta
lambda=m$lambda


margin(m$b0, m$beta, FHT$y, FHT$x, delta=2, loss = c("hsvm") )
KKT1(m$b0, m$beta, FHT$y, FHT$x, m$lambda,  ,thr=1e-03, delta=2, loss = c("hsvm"))
