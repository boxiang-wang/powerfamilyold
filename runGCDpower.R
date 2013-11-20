rm(list=ls(all=TRUE))

setwd("D:\\GitHub\\powerfamily")
require(Matrix)

source("GCDpower.R")
source("utilities.R")
source("p.GCDpower.R")



dyn.load("powerfamilyNET.dll")
dyn.load("hsvmlassoNET.dll")

load("FHT.rda")


# LASSO
m1 <- gcdnetpower(x=FHT$x,y=FHT$y,
                   #lambda=c(0.5,0.1),
                  lambda2=0, delta=0.01,method="hhsvm")
plot(m1, color=T)

m2 <- gcdnetpower(x=FHT$x,y=FHT$y,lambda2=0,qv=2,method="power")
plot(m2, color=T)

m3 <- gcdnetpower(x=FHT$x,y=FHT$y,
                  lambda2=0,qv=100,method="power")
plot(m3, color=T)

m4 <- gcdnetpower(x=FHT$x,y=FHT$y,        
                  lambda2=0,qv=0.01,method="power")
plot(m4, color=T)

print(predict.GCDpower(m1,type="class",newx=FHT$x[2:5,]))
print(predict.GCDpower(m2,type="class",newx=FHT$x[2:5,]))

