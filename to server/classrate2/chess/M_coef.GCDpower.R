######################################################################
## These functions are modifications from the
## gcdnet package:
## Yi Yang, Hui Zou, (2013).
## An Efficient Algorithm for Computing HHSVM and Its Generalization, 
## Journal of Computational and Graphical Statistics.
## http://users.stat.umn.edu/~zouxx019/Papers/gcdnet.pdf
## or http://users.stat.umn.edu/~yiyang/resources/papers/JCGS_gcdnet.pdf


coef.hsvmlassoNET <- function(object, s = NULL, type = c("coefficients", 
                                                     "nonzero"), ...) {
  type <- match.arg(type)
  b0 <- t(as.matrix(object$b0))
  rownames(b0) <- "(Intercept)"
  nbeta <- rbind2(b0, object$beta)
  if (!is.null(s)) {
    vnames <- dimnames(nbeta)[[1]]
    dimnames(nbeta) <- list(NULL, NULL)
    lambda <- object$lambda
    lamlist <- lambda.interp(lambda, s)
    nbeta <- nbeta[, lamlist$left, drop = FALSE] * lamlist$frac + 
      nbeta[, lamlist$right, drop = FALSE] * (1 - lamlist$frac)
    dimnames(nbeta) <- list(vnames, paste(seq(along = s)))
  }
  if (type == "coefficients") 
    return(nbeta)
  if (type == "nonzero") 
    return(nonzero(nbeta[-1, , drop = FALSE], bystep = TRUE))
} 
coef.powerfamily <- function(object, s = NULL, type = c("coefficients", 
                                                     "nonzero"), ...) {
  type <- match.arg(type)
  b0 <- t(as.matrix(object$b0))
  rownames(b0) <- "(Intercept)"
  nbeta <- rbind2(b0, object$beta)
  if (!is.null(s)) {
    vnames <- dimnames(nbeta)[[1]]
    dimnames(nbeta) <- list(NULL, NULL)
    lambda <- object$lambda
    lamlist <- lambda.interp(lambda, s)
    nbeta <- nbeta[, lamlist$left, drop = FALSE] * lamlist$frac + 
      nbeta[, lamlist$right, drop = FALSE] * (1 - lamlist$frac)
    dimnames(nbeta) <- list(vnames, paste(seq(along = s)))
  }
  if (type == "coefficients") 
    return(nbeta)
  if (type == "nonzero") 
    return(nonzero(nbeta[-1, , drop = FALSE], bystep = TRUE))
} 

coef.cv.GCDpower <- function(object, s = c("lambda.1se", 
                                         "lambda.min"), ...) {
  if (is.numeric(s)) 
    lambda <- s else if (is.character(s)) {
      s <- match.arg(s)
      lambda <- object[[s]]
    } else stop("Invalid form for s")
  coef(object$GCDpower.fit, s = lambda, ...)
} 


coef.GCDpower <- function(object, s = NULL, type = c("coefficients", 
                                                     "nonzero"), ...) NextMethod("coef") 
