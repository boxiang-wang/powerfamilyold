###############################################################################
# Modified old function
rm(list=ls(all=T))

load("D:/GitHub/powerfamily/data/SPECTF.rda")
require(gcdnet)
x = SPECTF.train[,-1]
y = c(-1, 1)[as.factor(SPECTF.train[, 1])]


setwd("D:\\GitHub\\powerfamily\\DrSVM")


# Initialize the parameters in the function
nfolds = 5
delta = 2
pred.loss = "misclass"

N <- nrow(x)
y <- drop(y)

# Fit a model first, in order to get lambda sequence
GCDpower.object <- gcdnet(x, y, lambda = NULL, delta = delta)
lambda <- GCDpower.object$lambda

# record the number of non-zero coefficients for each lambda
nz <- sapply(coef(GCDpower.object, type = "nonzero"), length)

# Start to initialize folder
set.seed(124)
foldid <- sample(rep(seq(nfolds), length = N)) 
outlist <- as.list(seq(nfolds))

# Within each folder, fit a model using training set
for (i in seq(nfolds)) {
  which <- foldid == i
  y_sub <- y[!which]
  outlist[[i]] <- gcdnet(x = x[!which, , drop = FALSE], 
                        y = y_sub, lambda = lambda, delta = delta, method="hhsvm")
}

# Compute the misclassification rate
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

# nfold rows and l column
# For each lambda, record misclassification rate for each folder
cvraw = (y != ifelse(predmat > 0, 1, -1))
outmat <- matrix(NA, nfolds, ncol(cvraw))
good <- matrix(0, nfolds, ncol(cvraw))
cvraw[is.infinite(cvraw)] <- NA
for (i in seq(nfolds)) {
  cvrawi <- cvraw[foldid == i, ]
  outmat[i, ] <- apply(cvrawi, 2, mean, na.rm = TRUE)
  good[i, seq(nlams[i])] <- 1
}
cvraw = outmat

# For each folder, record the length of lambda
N <- apply(good, 2, sum)
cvm <- apply(cvraw, 2, mean, na.rm = TRUE)
cvsd <- sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, mean, na.rm = TRUE)/(N - 1))
out <- list(lambda = lambda, cvm = cvm, cvsd = cvsd, 
            cvupper = cvm + cvsd, cvlo = cvm - cvsd, 
            nzero = nz, name = "Misclassification Error", 
            GCDpower.fit = GCDpower.object)
#lamin <- getmin(lambda, cvm, cvsd)

# Get the appropriate lambda
cvmin <- min(cvm)
idmin <- cvm <= cvmin
lambda.min <- max(lambda[idmin])
idmin <- match(lambda.min, lambda)
semin <- (cvm + cvsd)[idmin]
idmin <- cvm <= semin
lambda.1se <- max(lambda[idmin])
lamin <- list(lambda.min = lambda.min, lambda.1se = lambda.1se)
obj <- c(out, as.list(lamin))
class(obj) <- "cv.GCDpower"



#############################################################################
# Test new function cvs.GCDpower
rm(list=ls(all=T))

load("D:/GitHub/powerfamily/data/SPECTF.rda")
# require(gcdnet)
x = SPECTF.train[,-1]
y = c(-1, 1)[as.factor(SPECTF.train[, 1])]

setwd("D:\\GitHub\\powerfamily")
dyn.load("M_powerfamilyNET.dll")
source("M_GCDpower.R")
source("M_p.GCDpower.R")

setwd("D:\\GitHub\\powerfamily\\DrSVM")


cvs.GCDpower <- function(x, y, lambda = NULL, nfolds = 5, foldid, delta = 2, qv = 2, ...)
{
  
  # Initialize useful parameters
  pred.loss = "misclass"
  
  
  N <- nrow(x)
  y <- c(-1, 1)[as.factor(drop(y))]
  
  # Fit a model first, in order to get lambda sequence
  GCDpower.object <- gcdnetpower(x, y, lambda = lambda, delta = delta, qv = qv, 
                                ...)
  lambda <- GCDpower.object$lambda
  
  # record the number of non-zero coefficients for each lambda
  nz <- sapply(coef(GCDpower.object, type = "nonzero"), length)
  
  # Start to initialize folder
  foldid <- sample(rep(seq(nfolds), length = N)) 
  outlist <- as.list(seq(nfolds))             
  
  # Setting some parameters
  predmat <- matrix(NA, length(y), length(lambda)) # record prediction
  nlams <- double(nfolds)                   # record lambda length for each folder
  
  # Treat each folder as test set, and observations out of this folder 
  #   as training set. Fit models using training sets, and record the prediction
  #   on test sets.
  for (i in seq(nfolds)) {
    which <- foldid == i
    outlist[[i]] <- gcdnetpower(x = x[!which, , drop = FALSE], 
                           y = y[!which], lambda = lambda, delta = delta,
                                qv = qv, ...)
    preds <- predict(outlist[[i]], x[which, , drop = FALSE], type = "link")
    nlami <- length(outlist[[i]]$lambda)
    predmat[which, seq(nlami)] <- preds
    nlams[i] <- nlami
  }
  
  # For each lambda, record misclassification rate for each folder
  cvraw = (y != ifelse(predmat > 0, 1, -1))
  outmat <- matrix(NA, nfolds, ncol(cvraw))
  good <- matrix(0, nfolds, ncol(cvraw))
  cvraw[is.infinite(cvraw)] <- NA
  for (i in seq(nfolds)) {
    cvrawi <- cvraw[foldid == i, ]
    outmat[i, ] <- apply(cvrawi, 2, mean, na.rm = TRUE)
    good[i, seq(nlams[i])] <- 1
  }
  cvraw = outmat
  
  N <- apply(good, 2, sum)  # For each folder, record the length of lambda
  
  # Compute the mean and the se of the mean for each lambda's cv
  cvm <- apply(cvraw, 2, mean, na.rm = TRUE)
  cvsd <- sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, mean, na.rm = TRUE)/(N - 1))
  
  # Find the largest lambda that can achieve the minimum cv
  #   and the largest lambda than can achieve the minimum cv plus 1 se of cv
  cvmin <- min(cvm)
  idmin <- cvm <= cvmin
  lambda.min <- max(lambda[idmin])
  
  idmin <- match(lambda.min, lambda)
  semin <- (cvm + cvsd)[idmin]
  idmin <- cvm <= semin
  lambda.1se <- max(lambda[idmin])
  
  out <- list(lambda = lambda, cvm = cvm, cvsd = cvsd, 
              cvupper = cvm + cvsd, cvlo = cvm - cvsd, 
              nzero = nz, 
              lambda.min = lambda.min, lambda.1se = lambda.1se,
              name = "Misclassification Error", 
              GCDpower.fit = GCDpower.object)
  class(out) <- "cvs.GCDpower"
  out
}

set.seed(124)
m1 = cvs.GCDpower(x,y,method="hhsvm", delta=2)
m1$lambda.min

source("D:\\GitHub\\powerfamily\\M_cv.GCDpower.R")
set.seed(124)
m15 = cv.GCDpower(x,y,method="hhsvm", delta=2)
m15$lambda.min


set.seed(124)
require(gcdnet)
m2 = cv.gcdnet(x,y,method="hhsvm", delta=2, pred.loss="misclass")
m2$lambda.min



m1 = cv.gcdnet(x, y, delta = 2, method="hhsvm")
m2 = cv.GCDpower(x, y, delta = 2, method="hhsvm")
