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
source("M_cv.GCDpower.R")
# coefficients
source("M_coef.GCDpower.R")
# KKT
source("U_KKTcheckings.R")
# Source file of data generator
source("M_FHTgen.R")
source("DrSVM/DrSVM_Fix2.R")


# dyn.unload("M_powerfamilyNET.dll")
# shell("del M_powerfamilyNET.dll M_powerfamilyNET.o")
# shell("del M_powerfamilyintNET.dll M_powerfamilyintNET.o")
# shell("del M_powerfamilyhalfNET.dll M_powerfamilyhalfNET.o")
# shell("Rcmd SHLIB M_powerfamilyNET.f90 M_powerfamilyintNET.f90 M_powerfamilyhalfNET.f90 O_auxiliary.f90 -o M_powerfamilyNET.dll")


dyn.load("M_powerfamilyNET.dll")

# A small data
dat = NULL
dat$x = matrix(rnorm(80*10),nrow=80)
dat$y = sign(dat$x[,1]+dat$x[,2]+dat$x[,3]+dat$x[,4])
# g = svmL1L2(x,y,2)
# newx = matrix(rnorm(40*10),nrow=40)
# newy = sign(x[,1]+x[,2]+x[,3]+x[,4])
# pre = svmL1L2.predict(g,newx,newy)
#  Plotting:
# np = dim(x)
#	n = np[1]
#	p = np[2]
#	horizon1 <- apply(abs(g$beta), 1, sum)
#	matplot(horizon1, cbind(g$beta,pre$error), type="n", xlab="|beta|", ylab="beta")
#	for(i in 1:p){
#		lines(horizon1, g$beta[,i], col=i+1, type="l")
#	}
#	lines(horizon1,pre$error,col=p+2,type="l")

# A large data
load("data/colon.rda")
dat = NULL
dat$x = as.matrix(colon.x)
dat$y = c(-1,1)[as.factor(colon.y)]

# A large data
load("data/leuk.rda")
dat = NULL
dat$x = as.matrix(x)
dat$y = c(-1,1)[as.factor(y)]

#########################################################################
# comparison between power family and huber (gcdnet)

qv = 5
start1 = Sys.time()
m1 = gcdnetpower(x=dat$x, y=dat$y,
                lambda2=1, qv=qv, method="power",eps=1e-8, maxit=3e7, standardize=F)
stop1 = Sys.time()
difftime(stop1, start1, units="secs")

plot(m1)
title(paste("q=", qv, sep=""))
pre1 = predict(m1, newx = dat$x, type = "class", s=NULL)
error1 = mean(dat$y != pre1) 


require(gcdnet)
dd = 0.01
start1 = Sys.time()
m2 = gcdnet(x=dat$x, y=dat$y,
                 lambda2=1, delta=dd, method="hhsvm",eps=1e-8, maxit=3e7, standardize=F)
stop1 = Sys.time()
difftime(stop1, start1, units="secs")

plot(m2)
title(paste("d=", dd, sep=""))
pre2 = predict(m2, newx = dat$x, type = "class", s=NULL)
error2 = mean(dat$y != pre2) 
sum(abs(pre1-pre2))/2
mean(abs(pre1-pre2))/2

#########################################################################
# comparison between huber (gcdnet) and DrSVM
np = dim(dat$x)
n = np[1]
p = np[2]
delta = 0.01
lambda2D = 1*n*delta

source("DrSVM/drhsvm.r")
system.time(
md <- DrHSVM(x=dat$x,dat$y,lambda2D,delta=delta,eps=1e-10,#max.steps=1e5,
            type="lasso",scale=F, trace=T)
)
pre1 <- sign(DrHSVM.predict(md,dat$x,dat$y)$fit)
error1 = mean(dat$y != pre1) 

s <- apply(abs(md$beta), 1, sum)

#par(mfrow=c(1,1))

matplot(s, cbind(md$beta), type="l", cex.lab = 1.5, lty=1, col = gray.colors(12, start = 0.05, 
                                                                             end = 0.7, gamma = 2.2),
        xlab=expression(paste("||",beta,"||",scriptscriptstyle(1))), ylab=expression(beta),main="HHSVM")
 for(i in 1:p)
   lines(s, md$beta[,i], col=i+1, lty=1, pch =10)

system.time(
m12 <- gcdnet(x=dat$x, y=dat$y, lambda = md$lambda/n/delta,
            lambda2=1, delta=delta, method="hhsvm",eps=1e-10, maxit=1e5, standardize=F)
)
plot(m12, color=F)
for(i in 1:p)
  lines(s, m12$beta[i,], col=i+1, lty=1, pch =10)



pre12 = predict(m12, newx = dat$x, type = "class", s=NULL)
error12 = mean(dat$y != pre12) 
sum(abs(pre1-pre12))/2
mean(abs(pre1-pre12))/2
max(abs(t(md$beta) - m12$beta))

#########################################################################
# comparison between huber (gcdnet) and DrSVM_Fix2
np = dim(dat$x)
n = np[1]
p = np[2]
delta = 0.01
lambda2D = 1
source("DrSVM/DrSVM_Fix2.R")
system.time(
  md2 <- DrSVM_Fix2(x=dat$x, dat$y, lambda2D, eps=1e-10, #max_steps=3e7,
             smallmove=1,scale=F)
)


pre1 <- sign(t(DrSVM_Fix2.predict(md2,dat$x,dat$y)$f))
error1 = mean(dat$y != pre1) 

s <- apply(abs(md2$beta), 1, sum)
mm = cbind(md2$beta)
matplot(s, mm, type="l", cex.lab = 1.5, 
        lty=1, col = gray.colors(12, start = 0.05,
        end = 0.7, gamma = 2.2), xlim=c(0,4.5), ylim=c(-0.5,0.5),
        xlab=expression(paste("||",beta,"||",scriptscriptstyle(1))), 
        ylab=expression(beta),main="HHSVM")
for(i in 1:p)
  lines(s, md2$beta[,i], col=i+1, lty=1, pch =10)


m0001 <- gcdnet(x=dat$x, y=dat$y, lambda = md2$lambda/n,
            lambda2=lambda2D/n, delta=0.0001, method="hhsvm",eps=1e-10, maxit=3e7, standardize=F)

plot(m0001, xlim=c(0,4.5), ylim=c(-0.5,0.5), main="delta=0.0001")

m001 = gcdnet(x=dat$x, y=dat$y, lambda = md2$lambda/n,
             lambda2=lambda2D/n, delta=0.001, method="hhsvm",eps=1e-10, maxit=3e7, standardize=F)
plot(m001, xlim=c(0,4.5), ylim=c(-0.5,0.5), main="delta=0.001")

m01 = gcdnet(x=dat$x, y=dat$y, lambda = md2$lambda/n,
            lambda2=lambda2D/n, delta=0.01, method="hhsvm",eps=1e-10, maxit=3e7, standardize=F)
plot(m01, xlim=c(0,4.5), ylim=c(-0.5,0.5), main="delta=0.01")

m1 <- gcdnet(x=dat$x, y=dat$y, lambda = md2$lambda/n,
             lambda2=lambda2D/n, delta=0.1, method="hhsvm",eps=1e-10, maxit=3e7, standardize=F)
plot(m1, xlim=c(0,4.5), ylim=c(-0.5,0.5), main="delta=0.1")

ad = apply(abs(md2$beta), 1, sum)
a0001 = apply(abs(m0001$beta), 2, sum)
a001 = apply(abs(m001$beta), 2, sum)
a01 = apply(abs(m01$beta), 2, sum)
a1 = apply(abs(m1$beta), 2, sum)
cbind(ad,a0001,a001,a01,a1)

pre2 = predict(m2, newx = dat$x, type = "class", s=NULL)
error2 = mean(dat$y != pre2) 
mean(abs(pre1-pre2))/2
dim(abs(pre1-pre2))

rownames(md2$beta)=NULL
rownames(m001$beta)=NULL
max(abs(t(md2$beta) - m0001$beta))
max(abs(t(md2$beta) - m001$beta))
max(abs(t(md2$beta) - m01$beta))
max(abs(t(md2$beta) - m1$beta))

mean(abs(m01$beta - m1$beta))
mean(abs(m001$beta - m1$beta))
mean(abs(m0001$beta - m1$beta))
mean(abs(m01$beta - m1$beta))


m2=m1
pre2 = predict(m2, newx = dat$x, type = "class", s=NULL)
(error2 = mean(dat$y != pre2) )

#########################################################################
# comparison between huber (gcdnet) and DrSVM_Fix2 PLOTS
np = dim(dat$x)
n = np[1]
p = np[2]
delta = 0.01
lambda2D = 1
source("DrSVM/DrSVM_Fix2.R")
system.time(
  md2 <- DrSVM_Fix2(x=dat$x, dat$y, lambda2D, eps=1e-10, #max_steps=3e7,
                    smallmove=1,scale=F)
)

m001 = gcdnet(x=dat$x, y=dat$y, lambda = md2$lambda/n,
              lambda2=lambda2D/n, delta=0.001, method="hhsvm",eps=1e-10, maxit=3e7, standardize=F)
m01 = gcdnet(x=dat$x, y=dat$y, lambda = md2$lambda/n,
             lambda2=lambda2D/n, delta=0.01, method="hhsvm",eps=1e-10, maxit=3e7, standardize=F)
m1 <- gcdnet(x=dat$x, y=dat$y, lambda = md2$lambda/n,
             lambda2=lambda2D/n, delta=0.1, method="hhsvm",eps=1e-10, maxit=3e7, standardize=F)
solnpathplot = function(betas, lambda1, xlab = "L1 Norm", ylab = "Coefficients"
                        ,xlim=c(0,4.5), ylim=c(-0.5,0.5), main="plot")
{
  if(length(lambda1) != nrow(betas) && length(lambda1) == ncol(betas))
  {
    betas = t(betas) 
  }
  
  index <- apply(abs(betas), 1, sum)
  index <- c(0, index)
  betas <- rbind(0, as.matrix(betas))
  matplot(index, betas, 
          lty = 1, type = "l", 
          pch = 500, col = gray.colors(12, 
                                       start = 0.05, end = 0.7, gamma = 2.2),
          xlim=xlim, ylim=ylim,
          xlab=xlab, ylab=ylab, main=main) 
}
solnpathplot(md2$beta, md2$lambda, main="SVM")
solnpathplot(m001$beta, m0001$lambda, main="delta=0.0001")
solnpathplot(m01$beta, m01$lambda, main="delta=0.01")

#########################################################################
# cross validation
index = sample(1:nrow(dat$x),as.integer(nrow(dat$x)/3),replace=F) 
test_x = dat$x[index,]
test_y = dat$y[index]
train_x = dat$x[-index,]
train_y = dat$y[-index]


m1.cv = cv.GCDpower(x=train_x, y=train_y, lambda = md$lambda/n/delta,
                 lambda2=1, qv=50, method="power",eps=1e-8, maxit=3e7, standardize=F)
m1.cv$lambda.min
pre1 = predict(m1.cv$GCDpower.fit, newx = test_x, s = m1.cv$lambda.min, type = "class")
error1 = mean(test_y != pre1)
plot(m1.cv)

m2.cv = cv.gcdnet(x=train_x, y=train_y, lambda = md$lambda/n/delta, pred.loss="misclass", 
            lambda2=1, delta=dd, method="hhsvm",eps=1e-8, maxit=3e7, standardize=F)
m2.cv$lambda.min
pre2 = predict(m2.cv$gcdnet.fit, newx = test_x, s = m2.cv$lambda.min, type = "class")
error2 = mean(test_y != pre2) 
plot(m2.cv)