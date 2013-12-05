setwd("D:\\GitHub\\powerfamily")

load("D_FHT.rda")
y = FHT$y
x = FHT$x

# data setup
y <- drop(y)
x <- as.matrix(x)
np <- dim(x)
nobs <- as.integer(np[1])
nvars <- as.integer(np[2])
vnames <- colnames(x)
if (is.null(vnames)) 
  vnames <- paste("V", seq(nvars), sep = "")
if (length(y) != nobs) 
  stop("x and y have different number of observations")

# from default parameter
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


## Source files of tool functions
source("O_plot.gcdnet.R")
source("O_utilities.R")
require(Matrix)



# About FORTRAN subroutines
dyn.unload("O_hsvmlassoNET.dll")
dyn.unload("O_sqsvmlassoNET.dll")
dyn.unload("M_powerfamilyNET.dll")
dyn.unload("M_powerfamilyintNET.dll")

## cmd
shell("del O_hsvmlassoNET.dll O_hsvmlassoNET.o")
shell("Rcmd SHLIB O_hsvmlassoNET.f90 O_auxiliary.f90 -o O_hsvmlassoNET.dll")


shell("del O_sqsvmlassoNET.dll O_sqsvmlassoNET.o")
shell("Rcmd SHLIB O_sqsvmlassoNET.f90 O_auxiliary.f90 -o O_sqsvmlassoNET.dll")

shell("del M_powerfamilyNET.dll M_powerfamilyNET.o")
shell("Rcmd SHLIB M_powerfamilyNET.f90 O_auxiliary.f90 -o M_powerfamilyNET.dll")

shell("del M_powerfamilyintNET.dll M_powerfamilyintNET.o")
shell("Rcmd SHLIB M_powerfamilyintNET.f90 O_auxiliary.f90 -o M_powerfamilyintNET.dll")

dyn.load("O_hsvmlassoNET.dll")
dyn.load("O_sqsvmlassoNET.dll")
dyn.load("M_powerfamilyNET.dll")
dyn.load("M_powerfamilyintNET.dll")



start1 = Sys.time()
#################################################################################
qv = 2
qv = as.double(qv)
fit1 <- .Fortran("powerfamilyNET", qv, lam2, nobs, nvars, 
                 as.double(x), as.double(y), jd, pf, pf2, dfmax, pmax, nlam, 
                 flmin, ulam, eps, isd, maxit, nalam = integer(1), b0 = double(nlam), 
                 beta = double(pmax * nlam), ibeta = integer(pmax), nbeta = integer(nlam), 
                 alam = double(nlam), npass = integer(1), jerr = integer(1))

#################################################################################
# output
fit1 <- getoutput(fit1, maxit, pmax, nvars, vnames)
fit1 <- c(fit1, list(npasses = fit1$npass, jerr = fit1$jerr))
class(fit1) <- c("powerfamilyNET")
if (is.null(lambda)) 
  fit1$lambda <- lamfix(fit1$lambda)
#fit$call <- this.call
#################################################################################
class(fit1) <- c("gcdnet", class(fit1))
#fit1

plot.gcdnet(fit1)

stop1 = Sys.time()
difftime(stop1, start1, units="secs")




start1 = Sys.time()
#################################################################################
qv = 2
qv = as.integer(qv)
fit1.int <- .Fortran("powerfamilyintNET", qv, lam2, nobs, nvars, 
                 as.double(x), as.double(y), jd, pf, pf2, dfmax, pmax, nlam, 
                 flmin, ulam, eps, isd, maxit, nalam = integer(1), b0 = double(nlam), 
                 beta = double(pmax * nlam), ibeta = integer(pmax), nbeta = integer(nlam), 
                 alam = double(nlam), npass = integer(1), jerr = integer(1))

#################################################################################
# output
fit1.int <- getoutput(fit1.int, maxit, pmax, nvars, vnames)
fit1.int <- c(fit1.int, list(npasses = fit1.int$npass, jerr = fit1.int$jerr))
class(fit1.int) <- c("powerfamilyNET")
if (is.null(lambda)) 
  fit1.int$lambda <- lamfix(fit1.int$lambda)
#fit1.int$call <- this.call
#################################################################################
class(fit1.int) <- c("gcdnet", class(fit1.int))
#fit1.int

plot.gcdnet(fit1.int)

stop1 = Sys.time()
difftime(stop1, start1, units="secs")





start0 = Sys.time()
#################################################################################
# call Fortran core
fit0 <- .Fortran("sqsvmlassoNET", lam2, nobs, nvars, 
                as.double(x), as.double(y), jd, pf, pf2, dfmax, pmax, nlam, 
                flmin, ulam, eps, isd, maxit, nalam = integer(1), b0 = double(nlam), 
                beta = double(pmax * nlam), ibeta = integer(pmax), nbeta = integer(nlam), 
                alam = double(nlam), npass = integer(1), jerr = integer(1))
#################################################################################
# output
fit0 <- getoutput(fit0, maxit, pmax, nvars, vnames)
fit0 <- c(fit0, list(npasses = fit0$npass, jerr = fit0$jerr))
class(fit0) <- c("sqsvmpath")
if (is.null(lambda)) 
  fit0$lambda <- lamfix(fit0$lambda)
#fit0$call <- this.call
#################################################################################
class(fit0) <- c("gcdnet", class(fit0))
#fit0

plot.gcdnet(fit0)
stop0 = Sys.time()
difftime(stop0, start0, units="secs")



start2 = Sys.time()
#################################################################################
# call Fortran core
delta <- as.double(delta)
fit2 <- .Fortran("hsvmlassoNET", delta, lam2, nobs, nvars, 
                as.double(x), as.double(y), jd, pf, pf2, dfmax, pmax, nlam, 
                flmin, ulam, eps, isd, maxit, nalam = integer(1), b0 = double(nlam), 
                beta = double(pmax * nlam), ibeta = integer(pmax), nbeta = integer(nlam), 
                alam = double(nlam), npass = integer(1), jerr = integer(1))
#################################################################################
# output
fit2 <- getoutput(fit2, maxit, pmax, nvars, vnames)
fit2 <- c(fit2, list(npasses = fit2$npass, jerr = fit2$jerr))
class(fit2) <- c("hsvmpath")
if (is.null(lambda)) 
  fit2$lambda <- lamfix(fit2$lambda)
#fit2$call <- this.call
#################################################################################
class(fit2) <- c("gcdnet", class(fit2))
#fit2

plot.gcdnet(fit2)
stop2 = Sys.time()
difftime(stop2, start2, units="secs")
