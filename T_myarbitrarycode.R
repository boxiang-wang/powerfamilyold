
setwd("D:\\GitHub\\powerfamily")

load("FHT.rda")
y = FHT$y
x = FHT$x
y <- drop(y)
x <- as.matrix(x)
np <- dim(x)
nobs <- as.integer(np[1])
nvars <- as.integer(np[2])
vnames <- colnames(x)

nlambda = 100
lambda.factor = ifelse(nobs < nvars, 0.01,1e-04)
lambda = NULL
lambda2 = 0.5
pf = rep(1, nvars)
pf2 = rep(1, nvars)
#exclude, 
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

source("plot.gcdnet.R")
source("utilities.R")
require(Matrix)

#dyn.load("auxiliary.dll")
dyn.load("hsvmlassoNET.dll")
dyn.load("sqsvmlassoNET.dll")
dyn.load("powerfamilyNET.dll")

dyn.unload("hsvmlassoNET.dll")
dyn.unload("sqsvmlassoNET.dll")
dyn.unload("powerfamilyNET.dll")

## cmd
del hsvmlassoNET.dll hsvmlassoNET.o
Rcmd SHLIB hsvmlassoNET.f90 auxiliary.f90 -o hsvmlassoNET.dll

del sqsvmlassoNET.dll sqsvmlassoNET.o
Rcmd SHLIB sqsvmlassoNET.f90 auxiliary.f90 -o sqsvmlassoNET.dll

del powerfamilyNET.dll powerfamilyNET.o
Rcmd SHLIB powerfamilyNET.f90 auxiliary.f90 -o powerfamilyNET.dll


delta=0.01
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
  fit$lambda <- lamfix(outlist$lambda)
#fit$call <- this.call
#################################################################################
class(fit) <- c("gcdnet", class(fit))
fit

plot.gcdnet(fit)




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





qv = 2
qv = as.double(qv)
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

plot.gcdnet(fit2)