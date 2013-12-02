rm(list=ls(all=TRUE))
setwd("D:\\GitHub\\powerfamily")

# Set the data set
load("D_FHT.rda")
y = FHT$y
x = FHT$x
y <- drop(y)
x <- as.matrix(x)

#x = matrix(c(1,0,-1,-1,1,0),3,2)
#y=c(-1,-1,1)


np <- dim(x)
nobs <- as.integer(np[1])
nvars <- as.integer(np[2])
vnames <- colnames(x)

nlambda = 100
lambda.factor = ifelse(nobs < nvars, 0.01,1e-04)
lambda = NULL
lambda2 = 0
pf = rep(1, nvars)
pf2 = rep(1, nvars)
dfmax = nvars + 1
pmax = min(dfmax * 1.2, nvars)
standardize = TRUE
eps = 1e-08
maxit = 1e+06
delta = 2
if (is.null(vnames)) 
  vnames <- paste("V", seq(nvars), sep = "")

#parameter setup
if (length(pf) != nvars) 
  stop("The size of L1 penalty factor must be same as the number of input variables")
if (length(pf2) != nvars) 
  stop("The size of L2 penalty factor must be same as the number of input variables")
if (lambda2 < 0) 
  stop("lambda2 must be non-negative")
maxit <- as.integer(maxit)
lam2 <- as.double(lambda2)
pf <- as.double(pf)
pf2 <- as.double(pf2)
isd <- as.integer(standardize)

isd = 0
isd = as.integer(isd)
eps <- as.double(eps)
dfmax <- as.integer(dfmax)
pmax <- as.integer(pmax)
jd <- as.integer(0)

#################################################################################
#lambda setup
nlam <- as.integer(nlambda)
if (is.null(lambda)) {
  if (lambda.factor >= 1) 
    stop("lambda.factor should be less than 1")
  flmin <- as.double(lambda.factor)
  ulam <- double(1)
} else {
  #flmin=1 if user define lambda
  flmin <- as.double(1)
  if (any(lambda < 0)) 
    stop("lambdas should be non-negative")
  ulam <- as.double(rev(sort(lambda)))
  nlam <- as.integer(length(lambda))
}

#################################################################################
#data setup
y <- as.factor(y)
y <- c(-1, 1)[as.numeric(y)]
if (!all(y %in% c(-1, 1))) 
  stop("y should be a factor with two levels")
if (delta < 0) 
  stop("delta must be non-negative")
delta <- as.double(delta)
delta=2
source("plot.gcdnet.R")
source("utilities.R")
require(Matrix)

#dyn.load("auxiliary.dll")

dyn.load("sqsvmlassoNET.dll")
dyn.load("powerfamilyNET.dll")
dyn.unload("powerfamilyNET.dll")

dyn.unload("hsvmlassoNET.dll")
#dyn.unload("hsvmlassoNET.dll")

# del hsvmlassoNETtestfunction.dll hsvmlassoNETtestfunction.o
# Rcmd SHLIB hsvmlassoNETtestfunction.f90 auxiliary.f90 -o hsvmlassoNETtestfunction.dll
# dyn.load("hsvmlassoNETtestfunction.dll")
# dyn.unload("hsvmlassoNETtestfunction.dll")

x = matrix(c(1,0,-1,-1,1,0),3,2)
y=c(-1,-1,1)
#################################################################################
# call Fortran core
fit <- .Fortran("hsvmlassoNET", delta, lam2, nobs, nvars, 
                as.double(x), as.double(y), jd, pf, pf2, dfmax, pmax, nlam, 
                flmin, ulam, eps, isd, maxit, nalam = integer(1), b0 = double(nlam), 
                beta = double(pmax * nlam), ibeta = integer(pmax), nbeta = integer(nlam), 
                alam = double(nlam), npass = integer(1), jerr = integer(1))
#################################################################################
# output
fit <- getoutput(fit, maxit, pmax, nvars, vnames)
fit <- c(fit, list(npasses = fit$npass, jerr = fit$jerr))
class(fit) <- c("hsvmpath")
if (is.null(lambda)) 
  fit$lambda <- lamfix(fit$lambda)
#fit$call <- this.call
#################################################################################
class(fit) <- c("gcdnet", class(fit))
fit

plot.gcdnet(fit)

print(predict(fit,type="class",newx=FHT$x[2:5,]))

b0 <- t(as.matrix(fit$b0))
rownames(b0) <- "(Intercept)"
nbeta <- rbind2(b0, fit$beta)

dim(nbeta)





#################################################################################
# call Fortran core
fit1 <- .Fortran("sqsvmlassoNET", lam2, nobs, nvars, 
                as.double(x), as.double(y), jd, pf, pf2, dfmax, pmax, nlam, 
                flmin, ulam, eps, isd, maxit, nalam = integer(1), b0 = double(nlam), 
                beta = double(pmax * nlam), ibeta = integer(pmax), nbeta = integer(nlam), 
                alam = double(nlam), npass = integer(1), jerr = integer(1))
#################################################################################
# output
fit1 <- getoutput(fit1, maxit, pmax, nvars, vnames)
fit1 <- c(fit1, list(npasses = fit1$npass, jerr = fit1$jerr))
class(fit1) <- c("hsvmpath")
if (is.null(lambda)) 
  fit1$lambda <- lamfix(outlist$lambda)
#fit$call <- this.call
#################################################################################
class(fit1) <- c("gcdnet", class(fit1))
fit1

plot.gcdnet(fit1)



delta=0.5
x = matrix(c(1,0,-1,-1,1,0),3,2)
y=c(-1,-1,1)
qv = 2
qv = as.double(qv)


#dyn.load("powerfamilyNETn.dll")
#dyn.unload("powerfamilyNETn.dll")

#del powerfamilyNETn.dll powerfamilyNETn.o
#Rcmd SHLIB powerfamilyNETn.f90 auxiliary.f90 -o powerfamilyNETn.dll


#dyn.load("powerfamilyNETn.dll")
qv=as.double(2)

fit2 <- .Fortran("powerfamilyNET", qv, lam2, nobs, nvars, 
                 as.double(x), as.double(y), jd, pf, pf2, dfmax, pmax, nlam, 
                 flmin, ulam, eps, isd, maxit, nalam = integer(1), b0 = double(nlam), 
                 beta = double(pmax * nlam), ibeta = integer(pmax), nbeta = integer(nlam), 
                 alam = double(nlam), npass = integer(1), jerr = integer(1))


#################################################################################
# output
fit2 <- getoutput(fit2, maxit, pmax, nvars, vnames)
fit2 <- c(fit2, list(npasses = fit2$npass, jerr = fit2$jerr))
if (is.null(lambda)) 
  fit2$lambda <- lamfix(outlist$lambda)

class(fit2) <- c("gcdnet", class(fit2))

fit2

print(predict(fit2,type="class",newx=FHT$x[2:5,]))




plot.gcdnet(fit2)

install.packages('gcdnet')
require('gcdnet')

x_log <- matrix(rnorm(100*10),100,10)
y_log <- sample(c(-1,1),100,replace=TRUE)
# LASSO
m <- gcdnet(x=x_log,y=y_log,lambda2=0,method="log")
plot(m)


#load("FHT.rda")
#y = FHT$y
#x = FHT$x
#y <- drop(y)
#x <- as.matrix(x)



data(promotergene)
install.packages("DWD")
library(DWD)
gene <- kdwd(type~.,data=spam,C=100,scaled=TRUE,cross=1)
gene@fitted

x=spam[,2:58]


rm(list=ls(all=TRUE))
setwd("D:\\GitHub\\powerfamily")
library(gcdnet)
data(FHT)
m = gcdnet(x=promotergene[,2:58],y=promotergene$Class,delta=2,method=c("hhsvm"))