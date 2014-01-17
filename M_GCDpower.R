######################################################################
## These functions are modifications from the
## gcdnet package:
## Yi Yang, Hui Zou, (2013).
## An Efficient Algorithm for Computing HHSVM and Its Generalization, 
## Journal of Computational and Graphical Statistics.
## http://users.stat.umn.edu/~zouxx019/Papers/gcdnet.pdf
## or http://users.stat.umn.edu/~yiyang/resources/papers/JCGS_gcdnet.pdf
powerfamilyintpath <- function(x, y, nlam, flmin, ulam, isd, 
                            eps, dfmax, pmax, jd, pf, pf2, maxit, lam2, qv, nobs, nvars, 
                            vnames) {
  #################################################################################
  #data setup
  y <- as.factor(y)
  y <- c(-1, 1)[as.numeric(y)]
  if (!all(y %in% c(-1, 1))) 
    stop("y should be a factor with two levels")
  if (qv < 0) 
    stop("delta must be non-negative")
  qv <- as.integer(qv)
  #################################################################################
  # call Fortran core
  fit <- .Fortran("powerfamilyintNET", qv, lam2, nobs, nvars, 
                  as.double(x), as.double(y), jd, pf, pf2, dfmax, pmax, nlam, 
                  flmin, ulam, eps, isd, maxit, nalam = integer(1), b0 = double(nlam), 
                  beta = double(pmax * nlam), ibeta = integer(pmax), nbeta = integer(nlam), 
                  alam = double(nlam), npass = integer(1), jerr = integer(1))
  #################################################################################
  # output
  outlist <- getoutput(fit, maxit, pmax, nvars, vnames)
  outlist <- c(outlist, list(npasses = fit$npass, jerr = fit$jerr))
  class(outlist) <- c("powerfamily")
  outlist
} 

powerfamilypath <- function(x, y, nlam, flmin, ulam, isd, 
                     eps, dfmax, pmax, jd, pf, pf2, maxit, lam2, qv, nobs, nvars, 
                     vnames) {
  #################################################################################
  #data setup
  y <- as.factor(y)
  y <- c(-1, 1)[as.numeric(y)]
  if (!all(y %in% c(-1, 1))) 
    stop("y should be a factor with two levels")
  if (qv < 0) 
    stop("delta must be non-negative")
  qv <- as.double(qv)
  #################################################################################
  # call Fortran core
  fit <- .Fortran("powerfamilyNET", qv, lam2, nobs, nvars, 
                   as.double(x), as.double(y), jd, pf, pf2, dfmax, pmax, nlam, 
                   flmin, ulam, eps, isd, maxit, nalam = integer(1), b0 = double(nlam), 
                   beta = double(pmax * nlam), ibeta = integer(pmax), nbeta = integer(nlam), 
                   alam = double(nlam), npass = integer(1), jerr = integer(1))
  #################################################################################
  # output
  outlist <- getoutput(fit, maxit, pmax, nvars, vnames)
  outlist <- c(outlist, list(npasses = fit$npass, jerr = fit$jerr))
  class(outlist) <- c("powerfamily")
  outlist
} 

powerfamilyhalfpath <- function(x, y, nlam, flmin, ulam, isd, 
                            eps, dfmax, pmax, jd, pf, pf2, maxit, lam2, qv, nobs, nvars, 
                            vnames) {
  #################################################################################
  #data setup
  y <- as.factor(y)
  y <- c(-1, 1)[as.numeric(y)]
  if (!all(y %in% c(-1, 1))) 
    stop("y should be a factor with two levels")
  if (qv < 0) 
    stop("delta must be non-negative")
  qv <- as.double(qv)
  #################################################################################
  # call Fortran core
  fit <- .Fortran("powerfamilyhalfNET", qv, lam2, nobs, nvars, 
                  as.double(x), as.double(y), jd, pf, pf2, dfmax, pmax, nlam, 
                  flmin, ulam, eps, isd, maxit, nalam = integer(1), b0 = double(nlam), 
                  beta = double(pmax * nlam), ibeta = integer(pmax), nbeta = integer(nlam), 
                  alam = double(nlam), npass = integer(1), jerr = integer(1))
  #################################################################################
  # output
  outlist <- getoutput(fit, maxit, pmax, nvars, vnames)
  outlist <- c(outlist, list(npasses = fit$npass, jerr = fit$jerr))
  class(outlist) <- c("powerfamily")
  outlist
} 

hsvmpath <- function(x, y, nlam, flmin, ulam, isd, 
                     eps, dfmax, pmax, jd, pf, pf2, maxit, lam2, delta, nobs, nvars, 
                     vnames) {
  #################################################################################
  #data setup
  y <- as.factor(y)
  y <- c(-1, 1)[as.numeric(y)]
  if (!all(y %in% c(-1, 1))) 
    stop("y should be a factor with two levels")
  if (delta < 0) 
    stop("delta must be non-negative")
  delta <- as.double(delta)
  #################################################################################
  # call Fortran core
  fit <- .Fortran("hsvmlassoNET", delta, lam2, nobs, nvars, 
                  as.double(x), as.double(y), jd, pf, pf2, dfmax, pmax, nlam, 
                  flmin, ulam, eps, isd, maxit, nalam = integer(1), b0 = double(nlam), 
                  beta = double(pmax * nlam), ibeta = integer(pmax), nbeta = integer(nlam), 
                  alam = double(nlam), npass = integer(1), jerr = integer(1))
  #################################################################################
  # output
  outlist <- getoutput(fit, maxit, pmax, nvars, vnames)
  outlist <- c(outlist, list(npasses = fit$npass, jerr = fit$jerr))
  class(outlist) <- c("hsvmlassoNET")
  outlist
} 


gcdnetpower <- function(x, y, nlambda = 100, method = c("power","hhsvm"), 
                        lambda.factor = ifelse(nobs < nvars, 0.01,1e-04), lambda = NULL, 
                        lambda2 = 0, pf = rep(1, nvars), pf2 = rep(1, nvars), exclude, 
                   dfmax = nvars + 1, pmax = min(dfmax * 1.2, nvars), standardize = TRUE, 
                   eps = 1e-08, maxit = 1e+06, delta = 2, qv = 2) {
  #################################################################################
  #data setup
  method <- match.arg(method)
  this.call <- match.call()
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
  #################################################################################
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
  if (!missing(exclude)) {
    jd <- match(exclude, seq(nvars), 0)
    if (!all(jd > 0)) 
      stop("Some excluded variables out of range")
    jd <- as.integer(c(length(jd), jd))
  } else jd <- as.integer(0)
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
  if ( (method == "power") && (abs(qv %% 1) < eps))
  {
    method = "powerint"
    qv = as.integer(qv)
  }
  if (method == "power")
  {
    if(abs(qv %% 1) < eps)
    {
      method = "powerint"
      qv = as.integer(qv)
    } 
    if (abs(qv - 0.5) < eps)
    {
      method = "powerhalf"
    }  
  }
  #################################################################################
  fit <- switch(method, 
                powerhalf = powerfamilyhalfpath(x, y, nlam, flmin, 
                                             ulam, isd, eps, dfmax, pmax, jd, pf, pf2, 
                                                 maxit, lam2, qv, nobs, nvars, vnames),
                powerint = powerfamilyintpath(x, y, nlam, flmin, 
                                  ulam, isd, eps, dfmax, pmax, jd, pf, pf2, maxit, 
                                  lam2, qv, nobs, nvars, vnames),
                power = powerfamilypath(x, y, nlam, flmin, 
                                 ulam, isd, eps, dfmax, pmax, jd, pf, pf2, maxit, 
                                        lam2, qv, nobs, nvars, vnames), 
                hhsvm = hsvmpath(x, y, nlam, flmin, 
                                 ulam, isd, eps, dfmax, pmax, jd, pf, pf2, maxit, 
                                 lam2, delta, nobs, nvars, vnames)
                )
  if(method == "powerint") method = "power"
  if (is.null(lambda)) 
    fit$lambda <- lamfix(fit$lambda)
  fit$call <- this.call
  #################################################################################
  class(fit) <- c("GCDpower", class(fit))
  fit
} 

source("O_utilities.R")