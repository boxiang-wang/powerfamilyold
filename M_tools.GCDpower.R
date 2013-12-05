######################################################################
## These functions are modifications from the
## gcdnet package:
## Yi Yang, Hui Zou, (2013).
## An Efficient Algorithm for Computing HHSVM and Its Generalization, 
## Journal of Computational and Graphical Statistics.
## http://users.stat.umn.edu/~zouxx019/Papers/gcdnet.pdf
## or http://users.stat.umn.edu/~yiyang/resources/papers/JCGS_gcdnet.pdf



powerfamily = function(u, qv)
{
  ifelse(u > (qv/(qv+1)), 
             1/(u^qv)*(qv^qv)/((qv+1)^(qv+1)), 1 - u )
}
  

dpowerfamily = function(u, qv)
{
  ifelse(u > (qv/(qv+1)), 
             -1/(u^(qv+1))*( (qv/(qv+1))^(qv+1) ), -1 )
}  
  

cv.GCDpower <- function(x, y, lambda = NULL, pred.loss = "misclass",
                        nfolds = 5, foldid, delta = 2, qv = 2, ...) {
  if (missing(pred.loss)) 
    pred.loss <- "misclass" else pred.loss <- match.arg(pred.loss)
  N <- nrow(x)
  ###Fit the model once to get dimensions etc of output
  y <- drop(y)
  GCDpower.object <- gcdnetpower(x, y, lambda = lambda, delta = delta, qv = qv, 
                          ...)
  lambda <- GCDpower.object$lambda
  # predict -> coef
  nz <- sapply(coef(GCDpower.object, type = "nonzero"), length)
  if (missing(foldid)) 
    foldid <- sample(rep(seq(nfolds), length = N)) else nfolds <- max(foldid)
  if (nfolds < 3) 
    stop("nfolds must be bigger than 3; nfolds=10 recommended")
  outlist <- as.list(seq(nfolds))
  ###Now fit the nfold models and store them
  for (i in seq(nfolds)) {
    which <- foldid == i
    y_sub <- y[!which]
    outlist[[i]] <- gcdnetpower(x = x[!which, , drop = FALSE], 
                           y = y_sub, lambda = lambda, delta = delta, qv = qv,...)
  }
  ###What to do depends on the pred.loss and the model fit
  fun <- paste("cv", class(GCDpower.object)[[2]], sep = ".")
  cvstuff <- do.call(fun, list(outlist, lambda, x, y, foldid, 
                               pred.loss, delta, qv))
  cvm <- cvstuff$cvm
  cvsd <- cvstuff$cvsd
  cvname <- cvstuff$name
  out <- list(lambda = lambda, cvm = cvm, cvsd = cvsd, cvupper = cvm + 
                cvsd, cvlo = cvm - cvsd, nzero = nz, name = cvname, GCDpower.fit = GCDpower.object)
  lamin <- getmin(lambda, cvm, cvsd)
  obj <- c(out, as.list(lamin))
  class(obj) <- fun
  obj
} 

cv.hsvmlassoNET <- function(outlist, lambda, x, y, foldid, 
                            pred.loss="misclass", delta, qv) {
  typenames <- c(misclass = "Misclassification Error", loss = "Margin Based Loss")
  ###Turn y into c(0,1)
  y <- as.factor(y)
  y <- c(-1, 1)[as.numeric(y)]
  nfolds <- max(foldid)
  predmat <- matrix(NA, length(y), length(lambda))
  nlams <- double(nfolds)
  for (i in seq(nfolds)) {
    which <- foldid == i
    fitobj <- outlist[[i]]
    preds <- predict(fitobj, x[which, , drop = FALSE], type = "link")
    nlami <- length(outlist[[i]]$lambda)
    predmat[which, seq(nlami)] <- preds
    nlams[i] <- nlami
  }
  cvraw <- switch(pred.loss, loss = 2 * hubercls(y * predmat, 
                                                 delta), misclass = (y != ifelse(predmat > 0, 1, -1)))
  cvob <- cvcompute(cvraw, foldid, nlams)
  cvraw <- cvob$cvraw
  N <- cvob$N
  cvm <- apply(cvraw, 2, mean, na.rm = TRUE)
  cvsd <- sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, mean, na.rm = TRUE)/(N - 
                                                                           1))
  list(cvm = cvm, cvsd = cvsd, name = typenames[pred.loss])
} 

cv.powerfamily <- function(outlist, lambda, x, y, foldid, 
                        pred.loss="misclass", delta, qv) {
  typenames <- c(misclass = "Misclassification Error", loss = "Margin Based Loss")
  ###Turn y into c(0,1)
  y <- as.factor(y)
  y <- c(-1, 1)[as.numeric(y)]
  nfolds <- max(foldid)
  predmat <- matrix(NA, length(y), length(lambda))
  nlams <- double(nfolds)
  for (i in seq(nfolds)) {
    which <- foldid == i
    fitobj <- outlist[[i]]
    preds <- predict(fitobj, x[which, , drop = FALSE], type = "link")
    nlami <- length(outlist[[i]]$lambda)
    predmat[which, seq(nlami)] <- preds
    nlams[i] <- nlami
  }
  cvraw <- switch(pred.loss, misclass = (y != ifelse(predmat > 0, 1, -1)))
  cvob <- cvcompute(cvraw, foldid, nlams)
  cvraw <- cvob$cvraw
  N <- cvob$N
  cvm <- apply(cvraw, 2, mean, na.rm = TRUE)
  cvsd <- sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, mean, na.rm = TRUE)/(N - 
                                                                           1))
  list(cvm = cvm, cvsd = cvsd, name = typenames[pred.loss])
} 


predict.cv.GCDpower <- function(object, newx, s = c("lambda.1se", 
                                                  "lambda.min"), ...) {
  if (is.numeric(s)) 
    lambda <- s else if (is.character(s)) {
      s <- match.arg(s)
      lambda <- object[[s]]
    } else stop("Invalid form for s")
  predict(object$GCDpower.fit, newx, s = lambda, ...)
} 
