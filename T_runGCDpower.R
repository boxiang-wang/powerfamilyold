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

# Two FORTRAN subroutines.
dyn.load("M_powerfamilyNET.dll")
dyn.load("O_hsvmlassoNET.dll")

# Data files
load("D_FHT.rda")


# LASSO
m1 <- gcdnetpower(x=FHT$x,y=FHT$y,
                   #lambda=c(0.5,0.1),
                  lambda2=0, delta=0.01,method="hhsvm")
pdf("lsolnpth_power.pdf",9,6)
plot(m1, color=T)
dev.off()

m_100 <- gcdnetpower(x=FHT$x,y=FHT$y,lambda2=0,qv=100,method="power")
pdf("lsolnpth_power100.pdf",9,6)
plot(m_100, color=T, main="q=100",xvar="lambda")
dev.off()

m_10 <- gcdnetpower(x=FHT$x,y=FHT$y,lambda2=0,qv=10,method="power")
pdf("lsolnpth_power10.pdf",9,6)
plot(m_10, color=T, main="q=10",xvar="lambda")
dev.off()

m_5 <- gcdnetpower(x=FHT$x,y=FHT$y,lambda2=0,qv=5,method="power")
pdf("lsolnpth_power5.pdf",9,6)
plot(m_5, color=T, main="q=5",xvar="lambda")
dev.off()

m_2 <- gcdnetpower(x=FHT$x,y=FHT$y,lambda2=0,qv=2,method="power")
pdf("lsolnpth_power2.pdf",9,6)
plot(m_2, color=T, main="q=2",xvar="lambda")
dev.off()

m_1 <- gcdnetpower(x=FHT$x,y=FHT$y,lambda2=0,qv=1,method="power")
pdf("lsolnpth_power1.pdf",9,6)
plot(m_1, color=T, main="q=1",xvar="lambda")
dev.off()

m_1_l_05 <- gcdnetpower(x=FHT$x,y=FHT$y,lambda2=0.5,qv=1,method="power")
pdf("solnpth_power1_l05.pdf",9,6)
plot(m_1_l_05 , color=T, main="q=1, lambda=0.5")
dev.off()

m_1_l_1 <- gcdnetpower(x=FHT$x,y=FHT$y,lambda2=1,qv=1,method="power")
pdf("solnpth_power1_l1.pdf",9,6)
plot(m_1_l_1 , color=T, main="q=1, lambda=1")
dev.off()


m_05 <- gcdnetpower(x=FHT$x,y=FHT$y,lambda2=0,qv=0.5,method="power")
pdf("lsolnpth_power05.pdf",9,6)
plot(m_05, color=T, main="q=0.5",xvar="lambda")
dev.off()


m_01 <- gcdnetpower(x=FHT$x,y=FHT$y,lambda2=0,qv=0.1,method="power")
pdf("lsolnpth_power01.pdf",9,6)
plot(m_01, color=T, main="q=0.1",xvar="lambda")
dev.off()


m_001 <- gcdnetpower(x=FHT$x,y=FHT$y,lambda2=0,qv=0.01,method="power")
pdf("lsolnpth_power001.pdf",9,6)
plot(m_001, color=T, main="q=0.01",xvar="lambda")
dev.off()








m_100 <- gcdnetpower(x=FHT$x,y=FHT$y,lambda2=0,delta=100,method="hhsvm")
pdf("hsolnpth_power100.pdf",9,6)
plot(m_100, color=T, main="delta=100")
dev.off()

m_10 <- gcdnetpower(x=FHT$x,y=FHT$y,lambda2=0,delta=10,method="hhsvm")
pdf("hsolnpth_power10.pdf",9,6)
plot(m_10, color=T, main="delta=10")
dev.off()

m_5 <- gcdnetpower(x=FHT$x,y=FHT$y,lambda2=0,delta=5,method="hhsvm")
pdf("hsolnpth_power5.pdf",9,6)
plot(m_5, color=T, main="delta=5")
dev.off()

m_2 <- gcdnetpower(x=FHT$x,y=FHT$y,lambda2=0,delta=2,method="hhsvm")
pdf("hsolnpth_power2.pdf",9,6)
plot(m_2, color=T, main="delta=2")
dev.off()

m_1 <- gcdnetpower(x=FHT$x,y=FHT$y,lambda2=0,delta=1,method="hhsvm")
pdf("hsolnpth_power1.pdf",9,6)
plot(m_1, color=T, main="delta=1")
dev.off()

m_05 <- gcdnetpower(x=FHT$x,y=FHT$y,lambda2=0,delta=0.5,method="hhsvm")
pdf("hsolnpth_power05.pdf",9,6)
plot(m_05, color=T, main="delta=0.5")
dev.off()


m_01 <- gcdnetpower(x=FHT$x,y=FHT$y,lambda2=0,delta=0.1,method="hhsvm")
pdf("hsolnpth_power01.pdf",9,6)
plot(m_01, color=T, main="delta=0.1")
dev.off()


m_001 <- gcdnetpower(x=FHT$x,y=FHT$y,lambda2=0,delta=0.01,method="hhsvm")
pdf("hsolnpth_power001.pdf",9,6)
plot(m_001, color=T, main="delta=0.01")
dev.off()




m3 <- gcdnetpower(x=FHT$x,y=FHT$y,
                  lambda2=0,qv=100,method="power")
plot(m3, color=T)

m4 <- gcdnetpower(x=FHT$x,y=FHT$y,        
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
m1 <- gcdnetpower(x=x,y=dat$y,
                  #lambda=c(0.5,0.1),
                  lambda2=0, delta=2,method="hhsvm")
plot(m1, color=T)
stop1=Sys.time()
difftime(stop1, start1, units="secs")

start1=Sys.time()
m2 <- gcdnetpower(x=x,y=dat$y,
                  #lambda=c(0.5,0.1),
                  lambda2=0, qv=1,method="power")
plot(m2, color=T)
stop1=Sys.time()
difftime(stop1, start1, units="secs")

start1=Sys.time()
m1 <- gcdnetpower(x=train.x,y=train.y,
                  #lambda=c(0.5,0.1),
                  lambda2=0, delta=2, method="hhsvm",nlambda=10)
stop1=Sys.time()
difftime(stop1, start1, units="secs")
plot(m1, color=T)


start1=Sys.time()
m2 <- gcdnetpower(x=train.x,y=train.y,
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
  m = gcdnetpower(x=train.x,y=train.y,
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



m = gcdnetpower(x=train.x,y=train.y,
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
  m = gcdnetpower(x=train.x,y=train.y,
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



m = gcdnetpower(x=train.x,y=train.y,
                #lambda=lambda.seq,
                lambda2=0, delta=qv.seq[i], qv=qv.seq[i], method="power",nlambda=100)

pdf('spam_q_100.pdf', 9, 6)
plot(m)
dev.off()

m00 = gcdnetpower(x=train.x,y=train.y,
                #lambda=lambda.seq,
                lambda2=0, delta=qv.seq[i], qv=1, method="power",nlambda=100)

pdf('spampower_q_1.pdf', 9, 6)
plot(m00, color=T, main="q=1")
dev.off()


m01 = gcdnetpower(x=train.x,y=train.y,
                #lambda=lambda.seq,
                lambda2=1, delta=qv.seq[i], qv=1, method="power",nlambda=100)

pdf('spampower_q_1_l_2.pdf', 9, 6)
plot(m01, color=T, main="q=1, lambda2=1")
dev.off()
coef(m)



m1 = gcdnetpower(x=train.x,y=train.y,
                #lambda=lambda.seq,
                lambda2=0, delta=qv.seq[i], qv=100, method="power",nlambda=100)

pdf('spampower_q_100.pdf', 9, 6)
plot(m1, color=T, main="q=100")
dev.off()



m2 = gcdnetpower(x=train.x,y=train.y,
                 #lambda=lambda.seq,
                 lambda2=0, delta=0.01, qv=1, method="hhsvm",nlambda=100)

pdf('spamhhsvm_d_001.pdf', 9, 6)
plot(m2, color=T, main="delta=0.01")
dev.off()

