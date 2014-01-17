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

# dyn.unload("M_powerfamilyNET.dll")
# shell("del M_powerfamilyNET.dll M_powerfamilyNET.o")
# shell("del M_powerfamilyintNET.dll M_powerfamilyintNET.o")
# shell("del M_powerfamilyhalfNET.dll M_powerfamilyhalfNET.o")
# shell("Rcmd SHLIB M_powerfamilyNET.f90 M_powerfamilyintNET.f90 M_powerfamilyhalfNET.f90 O_auxiliary.f90 -o M_powerfamilyNET.dll")

# dyn.load("M_powerfamilyNET.dll")

set.seed(1234)
FHT = FHTgen(n=80, p=95, rho=0.5)
dat = FHT

# timing 

start1 = Sys.time()
m = gcdnetpower(x=dat$x, y=dat$y,
                lambda2=0.01, qv=2, method="power",eps=1e-8, standardize=F)
stop1 = Sys.time()
difftime(stop1, start1, units="secs")

start1 = Sys.time()
m1 = gcdnetpower(x=dat$x, y=dat$y,
                lambda2=0.01, qv=1, method="power",eps=1e-8, standardize=F)
stop1 = Sys.time()
difftime(stop1, start1, units="secs")


system.time(m <- gcdnetpower(x=dat$x, y=dat$y,
                            lambda2=0.01, qv=2, method="power",eps=1e-8, standardize=F))[3]

KKT(m$b0, m$beta, dat$y, dat$x, m$lambda, lambda2=1.5, thr=1e-03, 
    qv=0.5, loss = c("power"))


start1 = Sys.time()
m = gcdnetpower(x=dat$x, y=dat$y,
                lambda2=1.5, delta=2/9, method="hhsvm",eps=1e-8, standardize=F)
stop1 = Sys.time()
difftime(stop1, start1, units="secs")


lambda2 = 1
qv = 1

# cv
cvs<-cv.GCDpower(dat$x, dat$y, eps=1e-8, qv=qv, delta=2,
                lambda2=lambda2, method="power",
                pred.loss="misclass", nfolds=5)

cvs$nzero
plot(cvs)
coef(cvs, type="nonzero")
class(cvs)

coef.cvs = coef(cvs, s="lambda.1se")[-1,]
length(coef.cvs[coef.cvs != 0])

# nonzero beta
m = gcdnetpower(x=dat$x, y=dat$y,
                lambda2=1.5, delta=2/9, method="power",eps=1e-8, standardize=F)
apply(coef(m), 2, function(x) length(x[x!=0]))

m$lambda
coef(m, s=m$lambda[1])

m = gcdnetpower(x=dat$x, y=dat$y,
                lambda2=1e-4, delta=2/9, method="power",eps=1e-8, standardize=F)
sum(coef(m, s=0.1)!=0)

m = gcdnetpower(x=dat$x, y=dat$y,
                lambda2=1, delta=2/9, method="power",eps=1e-8, standardize=F)
sum(coef(m, s=0.1)!=0)

