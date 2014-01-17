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
lambda2 = 1
pf = rep(1, nvars)
pf2 = rep(1, nvars)
#exclude, 
dfmax = nvars + 1
pmax = min(dfmax * 1.2, nvars)
standardize = FALSE
eps = 1e-08
maxit = 1e+06
delta = 2
qv=0.5

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
#dyn.unload("O_hsvmlassoNET.dll")
#dyn.unload("O_sqsvmlassoNET.dll")
dyn.unload("M_powerfamilyNET.dll")
dyn.unload("M_powerfamilyintNET.dll")
dyn.unload("M_powerfamilyhalfNET.dll")

## cmd
#shell("del O_hsvmlassoNET.dll O_hsvmlassoNET.o")
#shell("Rcmd SHLIB O_hsvmlassoNET.f90 O_auxiliary.f90 -o O_hsvmlassoNET.dll")


#shell("del O_sqsvmlassoNET.dll O_sqsvmlassoNET.o")
#shell("Rcmd SHLIB O_sqsvmlassoNET.f90 O_auxiliary.f90 -o O_sqsvmlassoNET.dll")

shell("del M_powerfamilyNET.dll M_powerfamilyNET.o")
#shell("del M_powerfamilyintNET.dll M_powerfamilyintNET.o")
#shell("del M_powerfamilyhalfNET M_powerfamilyhalfNET.o")

shell("Rcmd SHLIB M_powerfamilyNET.f90 M_powerfamilyintNET.f90 M_powerfamilyhalfNET.f90 O_auxiliary.f90 -o M_powerfamilyNET.dll")

#dyn.load("O_hsvmlassoNET.dll")
#dyn.load("O_sqsvmlassoNET.dll")
dyn.load("M_powerfamilyNET.dll")
#dyn.load("M_powerfamilyintNET.dll")
#dyn.load("M_powerfamilyhalfNET.dll")

start1 = Sys.time()
#################################################################################
qv = 0.5
qv = as.double(qv)
fit1 <- .Fortran("powerfamilyNET", qv, lam2, nobs, nvars, 
                 as.double(x), as.double(y), jd, pf, pf2, dfmax, pmax, nlam, 
                 flmin, ulam, eps, isd, maxit, nalam = integer(1), b0 = double(nlam), 
                 beta = double(pmax * nlam), ibeta = integer(pmax), nbeta = integer(nlam), 
                 alam = double(nlam), npass = integer(1), jerr = integer(1))

#fit1=fit
#################################################################################
# output
outlst <- getoutput(fit1, maxit, pmax, nvars, vnames)
outlst <- c(outlst, list(npasses = fit1$npass, jerr = fit1$jerr))
fit1=outlst
class(fit1) <- c("powerfamilyNET")
if (is.null(lambda)) 
  fit1$lambda <- lamfix(fit1$lambda)
#fit$call <- this.call
#################################################################################
class(fit1) <- c("gcdnet", class(fit1))
#summary(fit1)

plot.gcdnet(fit1)

stop1 = Sys.time()
difftime(stop1, start1, units="secs")
m = fit1


start1 = Sys.time()
#################################################################################
qv = 0.5
qv = as.double(qv)
fit1 <- .Fortran("powerfamilyhalfNET", qv, lam2, nobs, nvars, 
                 as.double(x), as.double(y), jd, pf, pf2, dfmax, pmax, nlam, 
                 flmin, ulam, eps, isd, maxit, nalam = integer(1), b0 = double(nlam), 
                 beta = double(pmax * nlam), ibeta = integer(pmax), nbeta = integer(nlam), 
                 alam = double(nlam), npass = integer(1), jerr = integer(1))

#fit1=fit
#################################################################################
# output
outlst <- getoutput(fit1, maxit, pmax, nvars, vnames)
outlst <- c(outlst, list(npasses = fit1$npass, jerr = fit1$jerr))
fit1=outlst
class(fit1) <- c("powerfamilyNET")
if (is.null(lambda)) 
  fit1$lambda <- lamfix(fit1$lambda)
#fit$call <- this.call
#################################################################################
class(fit1) <- c("gcdnet", class(fit1))
#summary(fit1)

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


####################################################################
# This program aims to find the best lambda_1
rm(list=ls(all=TRUE))


setwd("D:\\GitHub\\powerfamily")
require(Matrix)

# Source files with tool functions.
source("O_utilities.R")

# Main program
source("M_GCDpower.R")

# Prediction, plot
source("M_p.GCDpower.R")
# KKT checking, CV
source("M_tools.GCDpower.R")
# coefficients
source("M_coef.GCDpower.R")
# KKT
source("U_KKTcheckings.R")

dyn.load("M_powerfamilyNET.dll")

qv = 1
l2 = 0
# Data files
load("D_FHT.rda")
dat=FHT
x=dat$x
y=dat$y
N=nrow(x)
m = GCDpower(x=x, y=y,
                lambda2 = l2, qv = qv, method="power")
lambda = m$lambda
# Parameters of the function
nfolds = 5
foldid <- sample(rep(seq(nfolds), length = N))
outlist <- as.list(seq(nfolds))
###Now fit the nfold models and store them
for (i in seq(nfolds)) {
  which <- foldid == i
  y_train <- y[!which]
  x_train <- x[!which, , drop = FALSE]
  outlist[[i]] <- GCDpower(x = x_train, y = y_train, lambda = lambda,
                              lambda2 = l2, qv = qv, method="power")
}

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
cvraw <- (y != ifelse(predmat > 0, 1, -1))
mat= cvraw
nfolds <- max(foldid)
outmat <- matrix(NA, nfolds, ncol(mat))
good <- matrix(0, nfolds, ncol(mat))
for (i in seq(nfolds)) {
  mati <- mat[foldid == i, ]
  outmat[i, ] <- apply(mati, 2, mean, na.rm = TRUE)
  good[i, seq(nlams[i])] <- 1
}
N <- apply(good, 2, sum)
cvob = list(cvraw = outmat, N = N)
cvraw <- cvob$cvraw
N <- cvob$N
cvm <- apply(cvraw, 2, mean, na.rm = TRUE)
cvsd <- sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, mean, na.rm = TRUE)/(N - 1))

cvstuff = list(cvm = cvm, cvsd = cvsd, name = "Misclassification Error")

cvm <- cvstuff$cvm
cvsd <- cvstuff$cvsd
cvname <- cvstuff$name
out <- list(lambda = lambda, cvm = cvm, cvsd = cvsd, cvupper = cvm + 
              cvsd, cvlo = cvm - cvsd, nzero = nz, name = cvname, GCDpower.fit = GCDpower.object)
#lamin <- getmin(lambda, cvm, cvsd)
cvmin <- min(cvm)
idmin <- cvm <= cvmin
lambda.min <- max(lambda[idmin])
idmin <- match(lambda.min, lambda)
semin <- (cvm + cvsd)[idmin]
idmin <- cvm <= semin

lambda.1se <- max(lambda[idmin])
list(lambda.min = lambda.min, lambda.1se = lambda.1se)

obj <- c(out, as.list(lamin))
class(obj) <- fun
obj




#######################################
# Old code







# LASSO
m1 <- GCDpower(x=FHT$x,y=FHT$y,
                  #lambda=c(0.5,0.1),
                  lambda2=0, delta=0.01,method="hhsvm")
pdf("lsolnpth_power.pdf",9,6)
plot(m1, color=T)
dev.off()

m_100 <- GCDpower(x=FHT$x,y=FHT$y,lambda2=0,qv=100,method="power")
pdf("lsolnpth_power100.pdf",9,6)
plot(m_100, color=T, main="q=100",xvar="lambda")
dev.off()

m_10 <- GCDpower(x=FHT$x,y=FHT$y,lambda2=0,qv=10,method="power")
pdf("lsolnpth_power10.pdf",9,6)
plot(m_10, color=T, main="q=10",xvar="lambda")
dev.off()

m_5 <- GCDpower(x=FHT$x,y=FHT$y,lambda2=0,qv=5,method="power")
pdf("lsolnpth_power5.pdf",9,6)
plot(m_5, color=T, main="q=5",xvar="lambda")
dev.off()

m_2 <- GCDpower(x=FHT$x,y=FHT$y,lambda2=0,qv=2,method="power")
pdf("lsolnpth_power2.pdf",9,6)
plot(m_2, color=T, main="q=2",xvar="lambda")
dev.off()

m_1 <- GCDpower(x=FHT$x,y=FHT$y,lambda2=0,qv=1,method="power")
pdf("lsolnpth_power1.pdf",9,6)
plot(m_1, color=T, main="q=1",xvar="lambda")
dev.off()

m_1_l_05 <- GCDpower(x=FHT$x,y=FHT$y,lambda2=0.5,qv=1,method="power")
pdf("solnpth_power1_l05.pdf",9,6)
plot(m_1_l_05 , color=T, main="q=1, lambda=0.5")
dev.off()

m_1_l_1 <- GCDpower(x=FHT$x,y=FHT$y,lambda2=1,qv=1,method="power")
pdf("solnpth_power1_l1.pdf",9,6)
plot(m_1_l_1 , color=T, main="q=1, lambda=1")
dev.off()


m_05 <- GCDpower(x=FHT$x,y=FHT$y,lambda2=0,qv=0.5,method="power")
pdf("lsolnpth_power05.pdf",9,6)
plot(m_05, color=T, main="q=0.5",xvar="lambda")
dev.off()


m_01 <- GCDpower(x=FHT$x,y=FHT$y,lambda2=0,qv=0.1,method="power")
pdf("lsolnpth_power01.pdf",9,6)
plot(m_01, color=T, main="q=0.1",xvar="lambda")
dev.off()


m_001 <- GCDpower(x=FHT$x,y=FHT$y,lambda2=0,qv=0.01,method="power")
pdf("lsolnpth_power001.pdf",9,6)
plot(m_001, color=T, main="q=0.01",xvar="lambda")
dev.off()








m_100 <- GCDpower(x=FHT$x,y=FHT$y,lambda2=0,delta=100,method="hhsvm")
pdf("hsolnpth_power100.pdf",9,6)
plot(m_100, color=T, main="delta=100")
dev.off()

m_10 <- GCDpower(x=FHT$x,y=FHT$y,lambda2=0,delta=10,method="hhsvm")
pdf("hsolnpth_power10.pdf",9,6)
plot(m_10, color=T, main="delta=10")
dev.off()

m_5 <- GCDpower(x=FHT$x,y=FHT$y,lambda2=0,delta=5,method="hhsvm")
pdf("hsolnpth_power5.pdf",9,6)
plot(m_5, color=T, main="delta=5")
dev.off()

m_2 <- GCDpower(x=FHT$x,y=FHT$y,lambda2=0,delta=2,method="hhsvm")
pdf("hsolnpth_power2.pdf",9,6)
plot(m_2, color=T, main="delta=2")
dev.off()

m_1 <- GCDpower(x=FHT$x,y=FHT$y,lambda2=0,delta=1,method="hhsvm")
pdf("hsolnpth_power1.pdf",9,6)
plot(m_1, color=T, main="delta=1")
dev.off()

m_05 <- GCDpower(x=FHT$x,y=FHT$y,lambda2=0,delta=0.5,method="hhsvm")
pdf("hsolnpth_power05.pdf",9,6)
plot(m_05, color=T, main="delta=0.5")
dev.off()


m_01 <- GCDpower(x=FHT$x,y=FHT$y,lambda2=0,delta=0.1,method="hhsvm")
pdf("hsolnpth_power01.pdf",9,6)
plot(m_01, color=T, main="delta=0.1")
dev.off()


m_001 <- GCDpower(x=FHT$x,y=FHT$y,lambda2=0,delta=0.01,method="hhsvm")
pdf("hsolnpth_power001.pdf",9,6)
plot(m_001, color=T, main="delta=0.01")
dev.off()




m3 <- GCDpower(x=FHT$x,y=FHT$y,
                  lambda2=0,qv=100,method="power")
plot(m3, color=T)

m4 <- GCDpower(x=FHT$x,y=FHT$y,        
                  lambda2=0,qv=1,method="power")
plot(m4, color=T)

print(predict.GCDpower(m1,type="class",newx=FHT$x[2:5,]))
print(predict.GCDpower(m2,type="class",newx=FHT$x[2:5,]))

library(DWD)
data(spam)
nobs.spam=nrow(spam)
spam$y=ifelse(spam$type == "spam",1,-1)

index=sample.int(nobs.spam,300)
dat = spam[index,]

x=dat[,1:57]

start1=Sys.time()
m1 <- GCDpower(x=x,y=dat$y,
                  #lambda=c(0.5,0.1),
                  lambda2=0, delta=2,method="hhsvm")
plot(m1, color=T)
stop1=Sys.time()
difftime(stop1, start1, units="secs")

start1=Sys.time()
m2 <- GCDpower(x=x,y=dat$y,
                  #lambda=c(0.5,0.1),
                  lambda2=0, qv=1,method="power")
plot(m2, color=T)
stop1=Sys.time()
difftime(stop1, start1, units="secs")



start1=Sys.time()
m1 <- GCDpower(x=train.x,y=train.y,
                  #lambda=c(0.5,0.1),
                  lambda2=0.01, qv=2, method="power")
stop1=Sys.time()
difftime(stop1, start1, units="secs")
plot(m1, color=T)


start1=Sys.time()
m2 <- GCDpower(x=train.x,y=train.y,
                  #lambda=c(0.5,0.1),
                  lambda2=0, qv=1, method="power",nlambda=10)
plot(m2, color=T)
stop1=Sys.time()
difftime(stop1, start1, units="secs")

pred1 = print(predict.GCDpower(m1,type="class",newx=test.x))
pred2 = print(predict.GCDpower(m2,type="class",newx=test.x))

colSums(pred1 == test.y)
colSums(m1$beta != 0)

colSums(pred2 == test.y)
colSums(m2$beta != 0)


start1=Sys.time()
m1 <- cv.GCDpower(x=train.x,y=train.y,
                  #lambda=c(0.5,0.1),
                  lambda2=0, delta=1.5, qv=1, method="hhsvm",nlambda=10)
stop1=Sys.time()
difftime(stop1, start1, units="secs")

pred1=predict.cv.GCDpower(m1, s="lambda.min",newx=test.x)



qv.seq=c(0.1,0.5,1,1.5,2,3,5)
qv.length=length(qv.seq)
misclassrate = rep(NA, qv.length)
nonzerobate = rep(NA, qv.length)
start1=Sys.time()
for(i in 1:qv.length)
{
  
  print(i)
  m = cv.GCDpower(x=train.x,y=train.y,
                  lambda2=0, delta=qv.seq[i], qv=qv.seq[i], method="hhsvm",nlambda=10)
  pred=predict.cv.GCDpower(m, s="lambda.min",newx=test.x)
  misclassrate[i] = as.numeric(colSums(pred == test.y))
  nonzerobate[i] = m$nzero[match(m$lambda.min, m$lambda)]
} 
cbind(qv.seq, misclassrate, nonzerobate)
stop1=Sys.time()
difftime(stop1, start1, units="secs")





set.seed(123)
load("FHT.rda")
index = sample.int(50,50)
folderind=rep(1:5,each=10)
train.index = index[folderind!=5]
test.index = index[folderind==5]

train.index=test.index=index
train.x = FHT$x[train.index,]; train.y = FHT$y[train.index]
test.x = FHT$x[test.index,]; test.y = FHT$y[test.index]


index = sample.int(50,50)
folderind=rep(1:5,each=10)
train.index = index[folderind!=5]
test.index = index[folderind==5]
train.x = FHT$x[train.index,]; train.y = FHT$y[train.index]
test.x = FHT$x[test.index,]; test.y = FHT$y[test.index]


qv.seq=c(0.01,0.5,1,2,5,100)
qv.length=length(qv.seq)

lambda.seq=c(0.1,0.15,0.2,0.3,0.5)
lambda.length=length(lambda.seq)
l = matrix(NA, qv.length, lambda.length)
misclassrate = matrix(NA, qv.length, lambda.length)
nonzerobeta = matrix(NA, qv.length, lambda.length)
start1=Sys.time()
for(i in 1:qv.length)
{
  
  print(i)
  m = GCDpower(x=train.x,y=train.y,
                  #lambda=lambda.seq,
                  lambda2=0, delta=qv.seq[i], qv=qv.seq[i], method="hhsvm",nlambda=100)
  for(j in 1:lambda.length)
  {
    pred=predict(m, lambda.seq[j],newx=test.x)
    l[i,j]=j
    misclassrate[i,j] = as.numeric(colSums(pred != test.y))/length(test.y)
    nonzerobeta[i,j] = sum(coef(m, s=lambda.seq[j])[-1,]!=0)
  }
} 

a=cbind(c(NA,qv.seq),rbind(lambda.seq,misclassrate))
stop1=Sys.time()
difftime(stop1, start1, units="secs")


b=cbind(c(NA,qv.seq),rbind(lambda.seq,nonzerobeta))



m = GCDpower(x=train.x,y=train.y,
                #lambda=lambda.seq,
                lambda2=0, delta=0.01, qv=10, method="hhsvm")
plot(m, color=F)


library(DWD)

set.seed(123)
index=sample.int(nobs.spam,300)
dat = spam[index,]

index = sample.int(300,300)
folderind=rep(1:5,each=60)
train.index = index[folderind!=5]
test.index = index[folderind==5]
train.x = dat[train.index,1:57]; train.y = dat$y[train.index]
test.x = dat[test.index,1:57]; test.y = dat$y[test.index]


i=6
qv.seq=c(0.01,0.5,1,2,5,100)

qv.seq=sort(qv.seq,T)
qv.length=length(qv.seq)

lambda.seq=c(0.1,0.15,0.2,0.3,0.5)
lambda.length=length(lambda.seq)
l = matrix(NA, qv.length, lambda.length)
misclassrate = matrix(NA, qv.length, lambda.length)
nonzerobeta = matrix(NA, qv.length, lambda.length)
start1=Sys.time()
for(i in 1:qv.length)
{
  
  print(i)
  m = GCDpower(x=train.x,y=train.y,
                  #lambda=lambda.seq,
                  lambda2=0, delta=qv.seq[i], qv=qv.seq[i], method="power",nlambda=100)
  for(j in 1:lambda.length)
  {
    pred=predict(m, lambda.seq[j],newx=test.x)
    l[i,j]=j
    misclassrate[i,j] = as.numeric(colSums(pred != test.y))/length(test.y)
    nonzerobeta[i,j] = sum(coef(m, s=lambda.seq[j])[-1,]!=0)
  }
} 

a=cbind(c(NA,qv.seq),rbind(lambda.seq,misclassrate*60))
b=cbind(c(NA,qv.seq),rbind(lambda.seq,nonzerobeta))
stop1=Sys.time()
difftime(stop1, start1, units="secs")

powera=a
powerb=b
powerdtime = difftime(stop1, start1, units="secs")

hhsvma = a
hhsvmb = b
hhsvmdtime = difftime(stop1, start1, units="secs")



m = GCDpower(x=train.x,y=train.y,
                #lambda=lambda.seq,
                lambda2=0, delta=qv.seq[i], qv=qv.seq[i], method="power",nlambda=100)

pdf('spam_q_100.pdf', 9, 6)
plot(m)
dev.off()

m00 = GCDpower(x=train.x,y=train.y,
                  #lambda=lambda.seq,
                  lambda2=0, delta=qv.seq[i], qv=1, method="power",nlambda=100)

pdf('spampower_q_1.pdf', 9, 6)
plot(m00, color=T, main="q=1")
dev.off()


m01 = GCDpower(x=train.x,y=train.y,
                  #lambda=lambda.seq,
                  lambda2=1, delta=qv.seq[i], qv=1, method="power",nlambda=100)

pdf('spampower_q_1_l_2.pdf', 9, 6)
plot(m01, color=T, main="q=1, lambda2=1")
dev.off()
coef(m)



m1 = GCDpower(x=train.x,y=train.y,
                 #lambda=lambda.seq,
                 lambda2=0, delta=qv.seq[i], qv=100, method="power",nlambda=100)

pdf('spampower_q_100.pdf', 9, 6)
plot(m1, color=T, main="q=100")
dev.off()



m2 = GCDpower(x=train.x,y=train.y,
                 #lambda=lambda.seq,
                 lambda2=0, delta=0.01, qv=1, method="hhsvm",nlambda=100)

pdf('spamhhsvm_d_001.pdf', 9, 6)
plot(m2, color=T, main="delta=0.01")
dev.off()



###############################
rm(list=ls(all=TRUE))

setwd("D:\\GitHub\\powerfamily")
require(Matrix)

# Source files with tool functions.
source("O_utilities.R")

# Main program
source("M_GCDpower.R")

# Prediction, plot
source("M_p.GCDpower.R")
# KKT checking, CV
source("M_tools.GCDpower.R")
# coefficients
source("M_coef.GCDpower.R")
# KKT
source("U_KKTcheckings.R")
dyn.load("M_powerfamilyNET.dll")

# Cross Validation
load("D_FHT.rda")
y = FHT$y
x = FHT$x

lambda=NULL

pred.loss <- "misclass"
nfolds = 5
qv = 2

N <- nrow(x)
###Fit the model once to get dimensions etc of output
y <- drop(y)
GCDpower.object <- GCDpower(x, y, lambda = lambda, 
                               qv = qv, 
                               method="power",
                               lambda2=0)
lambda <- GCDpower.object$lambda
# predict -> coef
nz <- sapply(coef(GCDpower.object, type = "nonzero"), length)

for (i in seq(nfolds)) {
  which <- foldid == i
  y_sub <- y[!which]
  outlist[[i]] <- GCDpower(x = x[!which, , drop = FALSE], 
                              y = y_sub, lambda = lambda,qv = qv, 
                              method="power",lambda2=0)
}

fun <- paste("cv", class(GCDpower.object)[[2]], sep = ".")
cvstuff <- do.call(fun, list(outlist, lambda, x, y, foldid, 
                             +                              pred.loss, delta, qv))