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


gcdnetpower <- function(x, y, nlambda = 100, method = c("power","hhsvm","logit", "sqsvm", "ls"), 
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


############################
## Source utilities file to get the formal object.
######################################################################
## These functions are minor modifications or directly
#   copied from the
## glmnet package:
## Jerome Friedman, Trevor Hastie, Robert Tibshirani
#   (2010).
## Regularization Paths for Generalized Linear Models via
#   Coordinate Descent.
##        Journal of Statistical Software, 33(1), 1-22.
##        URL http://www.jstatsoft.org/v33/i01/.
## The reason they are copied here is because they are
#   internal functions
## and hence are not exported into the global environment.
## The original comments and header are preserved.


cvcompute <- function(mat, foldid, nlams) {
  ###Computes the weighted mean and SD within folds, and
  #   hence
  #   the se of the mean
  nfolds <- max(foldid)
  outmat <- matrix(NA, nfolds, ncol(mat))
  good <- matrix(0, nfolds, ncol(mat))
  mat[is.infinite(mat)] <- NA
  for (i in seq(nfolds)) {
    mati <- mat[foldid == i, ]
    outmat[i, ] <- apply(mati, 2, mean, na.rm = TRUE)
    good[i, seq(nlams[i])] <- 1
  }
  N <- apply(good, 2, sum)
  list(cvraw = outmat, N = N)
}



err <- function(n, maxit, pmax) {
  if (n == 0) 
    msg <- ""
  if (n > 0) {
    if (n < 7777) 
      msg <- "Memory allocation error"
    if (n == 7777) 
      msg <- "All used predictors have zero variance"
    if (n == 10000) 
      msg <- "All penalty factors are <= 0"
    n <- 1
    msg <- paste("in gcdnet fortran code -", msg)
  }
  if (n < 0) {
    if (n > -10000) 
      msg <- paste("Convergence for ", -n, "th lambda value not reached after maxit=", 
                   maxit, " iterations; solutions for larger lambdas returned", 
                   sep = "")
    if (n < -10000) 
      msg <- paste("Number of nonzero coefficients along the path exceeds pmax=", 
                   pmax, " at ", -n - 10000, "th lambda value; solutions for larger lambdas returned", 
                   sep = "")
    n <- -1
    msg <- paste("from gcdnet fortran code -", msg)
  }
  list(n = n, msg = msg)
}



error.bars <- function(x, upper, lower, width = 0.02, 
                       ...) {
  xlim <- range(x)
  barw <- diff(xlim) * width
  segments(x, upper, x, lower, ...)
  segments(x - barw, upper, x + barw, upper, ...)
  segments(x - barw, lower, x + barw, lower, ...)
  range(upper, lower)
}


getmin <- function(lambda, cvm, cvsd) {
  cvmin <- min(cvm)
  idmin <- cvm <= cvmin
  lambda.min <- max(lambda[idmin])
  idmin <- match(lambda.min, lambda)
  semin <- (cvm + cvsd)[idmin]
  idmin <- cvm <= semin
  # cat('\n\nidmin\n\n',idmin)
  # cat('\n\nlambda[idmin]\n\n',lambda[idmin])
  # cat('\n\nmax\n\n',max(lambda[idmin]))
  lambda.1se <- max(lambda[idmin])
  list(lambda.min = lambda.min, lambda.1se = lambda.1se)
}


getoutput <- function(fit, maxit, pmax, nvars, vnames) {
  nalam <- fit$nalam
  nbeta <- fit$nbeta[seq(nalam)]
  nbetamax <- max(nbeta)
  lam <- fit$alam[seq(nalam)]
  stepnames <- paste("s", seq(nalam) - 1, sep = "")
  errmsg <- err(fit$jerr, maxit, pmax)
  switch(paste(errmsg$n), `1` = stop(errmsg$msg, call. = FALSE), 
         `-1` = print(errmsg$msg, call. = FALSE))
  dd <- c(nvars, nalam)
  if (nbetamax > 0) {
    beta <- matrix(fit$beta[seq(pmax * nalam)], pmax, nalam)[seq(nbetamax), 
                                                             , drop = FALSE]
    df <- apply(abs(beta) > 0, 2, sum)
    ja <- fit$ibeta[seq(nbetamax)]
    oja <- order(ja)
    ja <- rep(ja[oja], nalam)
    ibeta <- cumsum(c(1, rep(nbetamax, nalam)))
    beta <- new("dgCMatrix", Dim = dd, Dimnames = list(vnames, 
                                                       stepnames), x = as.vector(beta[oja, ]), p = as.integer(ibeta - 
                                                                                                                1), i = as.integer(ja - 1))
  } else {
    beta <- zeromat(nvars, nalam, vnames, stepnames)
    df <- rep(0, nalam)
  }
  b0 <- fit$b0
  if (!is.null(b0)) {
    b0 <- b0[seq(nalam)]
    names(b0) <- stepnames
  }
  list(b0 = b0, beta = beta, df = df, dim = dd, lambda = lam)
}



lambda.interp <- function(lambda, s) {
  ### lambda is the index sequence that is produced by the
  #   model
  ### s is the new vector at which evaluations are required.
  ### the value is a vector of left and right indices, and a
  #   vector of fractions.
  ### the new values are interpolated bewteen the two using
  #   the
  #   fraction
  ### Note: lambda decreases. you take:
  ### sfrac*left+(1-sfrac*right)
  if (length(lambda) == 1) {
    nums <- length(s)
    left <- rep(1, nums)
    right <- left
    sfrac <- rep(1, nums)
  } else {
    s[s > max(lambda)] <- max(lambda)
    s[s < min(lambda)] <- min(lambda)
    k <- length(lambda)
    sfrac <- (lambda[1] - s)/(lambda[1] - lambda[k])
    lambda <- (lambda[1] - lambda)/(lambda[1] - lambda[k])
    coord <- approx(lambda, seq(lambda), sfrac)$y
    left <- floor(coord)
    right <- ceiling(coord)
    sfrac <- (sfrac - lambda[right])/(lambda[left] - lambda[right])
    sfrac[left == right] <- 1
  }
  list(left = left, right = right, frac = sfrac)
}


lamfix <- function(lam) {
  llam <- log(lam)
  lam[1] <- exp(2 * llam[2] - llam[3])
  lam
}


nonzero <- function(beta, bystep = FALSE) {
  ns <- ncol(beta)
  ##beta should be in 'dgCMatrix' format
  if (nrow(beta) == 1) {
    if (bystep) 
      apply(beta, 2, function(x) if (abs(x) > 0) 
        1 else NULL) else {
          if (any(abs(beta) > 0)) 
            1 else NULL
        }
  } else {
    beta <- t(beta)
    which <- diff(beta@p)
    which <- seq(which)[which > 0]
    if (bystep) {
      nzel <- function(x, which) if (any(x)) 
        which[x] else NULL
      beta <- abs(as.matrix(beta[, which])) > 0
      if (ns == 1) 
        apply(beta, 2, nzel, which) else apply(beta, 1, nzel, which)
    } else which
  }
}



zeromat <- function(nvars, nalam, vnames, stepnames) {
  ca <- rep(0, nalam)
  ia <- seq(nalam + 1)
  ja <- rep(1, nalam)
  dd <- c(nvars, nalam)
  new("dgCMatrix", Dim = dd, Dimnames = list(vnames, stepnames), 
      x = as.vector(ca), p = as.integer(ia - 1), i = as.integer(ja - 
                                                                  1))
} 
