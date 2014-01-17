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
            lambda2=0, delta=2,method="hhsvm",eps=1e-6, standardize=F)
plot(m, color=T)

b0=m$b0
beta=m$beta
lambda=m$lambda


margin(m$b0, m$beta, FHT$y, FHT$x, delta=2, loss = c("hsvm") )
KKT1(m$b0, m$beta, FHT$y, FHT$x, m$lambda, thr=1e-03, delta=2, loss = c("hsvm"))


#################################################################################
############ checking KKT conditions for power family ################
#################################################################################
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

KKT1 = function(b0, beta, y, x, lambda, lambda2, thr, 
                loss = c("hsvm", "power"), delta=2, qv=2) {
  loss = match.arg(loss)
  dl = margin(b0, beta, y, x, loss=loss, delta=delta, qv=qv)
  ctr = 0
  ccounts = 0
  for (l in 1:length(lambda)) {
    p = nrow(beta)
    for(j in 1:p)
    {
      if(beta[j,l]==0)
      {
        BB = abs(dl[j,l]) - lambda[l]
        ccounts = ccounts + 1
        if (BB > thr) 
        {
          cat("violate at b = 0", BB, "\n")
          ctr <- ctr + 1
        }
      } else{
        AA = dl[j,l] + lambda[l] * sign(beta[j,l]) + lambda2 * beta[j,l]
        ccounts = ccounts + 1
        if (abs(sum(AA)) >= thr)
        {
          cat("violate at b != 0", abs(sum(AA)), "\n")
          ctr <- ctr + 1
        }
        
      }
    }
  }
  cat("# of violations", ctr/length(lambda), "\n")
  cat("% of violations", ctr/ccounts*100, "%", "\n")
  return(ctr/length(lambda))
} 



#library(gcdnet)
load("D_FHT.rda")
y = FHT$y
x = FHT$x

#x = matrix(c(1,0,-1,-1,1,0),3,2)
#y = c(-1,-1,1)

# Source files with tool functions.
source("O_utilities.R")

# Main program
source("M_GCDpower.R")

# Two FORTRAN subroutines.
dyn.load("M_powerfamilyNET.dll")
dyn.load("O_hsvmlassoNET.dll")



m <- GCDpower(x=FHT$x,y=FHT$y,
            #lambda=c(0.1,0.01),
            lambda2=1, qv=2, method="power",eps=1e-10, standardize=F)
plot(m, color=T)

b0=m$b0
beta=m$beta
lambda=m$lambda

# KKT1 = function(b0, beta, y, x, lambda, lambda2, thr, delta, loss = c("hsvm"))
# margin(m$b0, m$beta, FHT$y, FHT$x, delta=2, loss = c("power") )
KKT(m$b0, m$beta, FHT$y, FHT$x, m$lambda, lambda2=1, thr=1e-03, qv=2, loss = c("power"))



#################################################################################
############ for the use of source ################
#################################################################################
rm(list=ls(all=TRUE))
setwd("D:\\GitHub\\powerfamily")

require(Matrix)
# Source files with tool functions.
source("O_utilities.R")
# Main program
source("M_GCDpower.R")
# Source file of KKT
source("U_KKTcheckings.R")
# FORTRAN subroutines.
dyn.load("M_powerfamilyNET.dll")
# Source file of data generator
source("M_FHTgen.R")

set.seed(1234)
FHT = FHTgen(n=5000, p=100, rho=0.5)


dat = FHT


m = GCDpower(x=dat$x, y=dat$y,
                 lambda2=1, qv=0.5, method="power",eps=1e-8, standardize=F)


KKT(m$b0, m$beta, dat$y, dat$x, m$lambda, lambda2=1, thr=1e-4, 
                  qv=0.5, loss = c("power"), print.out=F)

#################################################################################
############ construct KKT tables ################
#################################################################################

set.seed(1234)
FHT = FHTgen(n=100, p=5000, rho=0.8)
dat = FHT

start1 = Sys.time()
KKTtb(dat, lambda2=0, qv=0.5, nm="n100p5Kr08l20q05")
stop1 = Sys.time() 
difftime(stop1, start1, units="secs")

start1 = Sys.time()
KKTtb(dat, lambda2=1, qv=0.5, nm="n100p5Kr08l21q05")
stop1 = Sys.time()
difftime(stop1, start1, units="secs")

start1 = Sys.time()
KKTtb(dat, lambda2=0, qv=1, nm="n100p5Kr08l20q1")
stop1 = Sys.time()
difftime(stop1, start1, units="secs")

start1 = Sys.time()
KKTtb(dat, lambda2=1, qv=1, nm="n100p5Kr08l21q1")
stop1 = Sys.time()
difftime(stop1, start1, units="secs")

start1 = Sys.time()
KKTtb(dat, lambda2=0, qv=2, nm="n100p5Kr08l20q2")
stop1 = Sys.time()
difftime(stop1, start1, units="secs")

start1 = Sys.time()
KKTtb(dat, lambda2=1, qv=2, nm="n100p5Kr08l21q2")
stop1 = Sys.time()
difftime(stop1, start1, units="secs")

#shell("shutdown -s -f -t 1")

file.loc = "D:\\GitHub\\powerfamily\\Outputs\\KKTrda\\"
nm="n100p5Kr08l20q05"
load(paste(file.loc, nm, ".rda", sep=""))

p1=perct.tb
p2=perct.tb
xtable(cbind(p1,p2))


set.seed(1234)
FHT = FHTgen(n=100, p=5000, rho=0.8)
dat = FHT
m <- GCDpower(x=dat$x,y=dat$y,
                 #lambda=c(0.1,0.01),
                 lambda2=1, qv=2, method="power",eps=1e-4, standardize=F)
KKT(m$b0, m$beta, dat$y, dat$x, m$lambda, lambda2=1, thr=10^(-5), 
    qv=2, loss = c("power"), print.out=F)


