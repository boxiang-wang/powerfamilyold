rm(list=ls(all=TRUE))

setwd("D:\\GitHub\\powerfamily")

source("GCDpower.R")
source("utilities.R")
source("plot.gcdnet.R")
require(Matrix)

dyn.load("powerfamilyNET.dll")
dyn.load("hsvmlassoNET.dll")

x_log <- matrix(rnorm(100*10),100,10)
y_log <- sample(c(-1,1),100,replace=TRUE)


# LASSO
m1 <- gcdnetpower(x=x_log,y=y_log,lambda2=0, delta=0.01,method="hhsvm")
plot(m1)

m2 <- gcdnetpower(x=x_log,y=y_log,lambda2=0,method="power")
plot(m2)

m3 <- gcdnetpower(x=x_log,y=y_log,lambda2=0,qv=100,method="power")
plot(m3)


