setwd("D:\\GitHub\\powerfamily")
source("O_utilities.R")
dyn.load("M_powerfamilyNET.dll")
dyn.load("M_powerfamilyintNET.dll")


source("M_FHTgen.R")
source("U_KKTcheckings.R")
source("M_GCDpower.R")


setwd("D:\\GitHub\\powerfamily\\to server")
set.seed(1234)
FHT = FHTgen(n=100, p=5000, rho=0.8)
dat = FHT



setwd("D:\\GitHub\\powerfamily\\to server\\l0q1")
save(dat, file="dat.rda")
l2 = 0
qvalue = 0.5


m <- GCDpower(x=dat$x,y=dat$y,
                 lambda2=l2, qv=qvalue, method="power",
                 eps=1e-6, standardize=F)
b0_6 = m$b0
betas_6 = m$beta
lambda_6 = m$lambda



save(b0_6, file="b0_6.rda")
save(betas_6, file="beta_6.rda")
save(lambda_6, file="lambda_6.rda")


m <- GCDpower(x=dat$x,y=dat$y,
                 lambda2=l2, qv=qvalue, method="power",
                 eps=1e-7, standardize=F)
b0_7 = m$b0
betas_7 = m$beta
lambda_7 = m$lambda


save(b0_7, file="b0_7.rda")
save(betas_7, file="beta_7.rda")
save(lambda_7, file="lambda_7.rda")

m <- GCDpower(x=dat$x,y=dat$y,
                 lambda2=l2, qv=qvalue, method="power",
                 eps=1e-8, standardize=F)
b0_8 = m$b0
betas_8 = m$beta
lambda_8 = m$lambda


save(b0_8, file="b0_8.rda")
save(betas_8, file="beta_8.rda")
save(lambda_8, file="lambda_8.rda")



m <- GCDpower(x=dat$x,y=dat$y,
                 lambda2=l2, qv=qvalue, method="power",
                 eps=1e-9, standardize=F)
b0_9 = m$b0
betas_9 = m$beta
lambda_9 = m$lambda


save(b0_9, file="b0_9.rda")
save(betas_9, file="beta_9.rda")
save(lambda_9, file="lambda_9.rda")

m <- GCDpower(x=dat$x,y=dat$y,
                 lambda2=l2, qv=qvalue, method="power",
                 eps=1e-10, standardize=F)
b0_10 = m$b0
betas_10 = m$beta
lambda_10 = m$lambda


save(b0_10, file="b0_10.rda")
save(betas_10, file="beta_10.rda")
save(lambda_10, file="lambda_10.rda")



#############################
#############################
start1 = Sys.time()

l2 = 0
qvalue = 1



load("dat.rda")

load("b0_6.rda")
load("beta_6.rda")
dim(betas_6)
dim(betas_6)
load("lambda_6.rda")

source("U_KKTcheckings.R")

kk = 2:5
KKT.perc = rep(NA, 4)
for(i in 1:length(kk))
{
  KKT.perc[i] = KKT(b0_6, betas_6, dat$y, dat$x, lambda_6, lambda2=l2, thr=10^(-kk[i]), 
                    qv=qvalue, loss = c("power"), print.out=F)
  write.csv(KKT.perc[i], file=paste("thr6_", kk[i], ".csv", sep=""))
}

write.csv(KKT.perc, file="KKT6.csv")
KKT.perc


load("b0_7.rda")
load("beta_7.rda")
load("lambda_7.rda")

dim(betas_7)
dim(betas_7)

kk = 2:5
KKT.perc = rep(NA, 4)
for(i in 1:length(kk))
{
  KKT.perc[i] = KKT(b0_7, betas_7, dat$y, dat$x, lambda_7, lambda2=l2, thr=10^(-kk[i]), 
                    qv=qvalue, loss = c("power"), print.out=F)
  write.csv(KKT.perc[i], file=paste("thr7_", kk[i], ".csv", sep=""))
}

write.csv(KKT.perc, file="KKT7.csv")
KKT.perc


load("b0_8.rda")
load("beta_8.rda")
load("lambda_8.rda")

dim(betas_8)
dim(betas_8)

kk = 2:5
KKT.perc = rep(NA, 4)
for(i in 1:length(kk))
{
  KKT.perc[i] = KKT(b0_8, betas_8, dat$y, dat$x, lambda_8, lambda2=l2, thr=10^(-kk[i]), 
                    qv=qvalue, loss = c("power"), print.out=F)
  write.csv(KKT.perc[i], file=paste("thr8_", kk[i], ".csv", sep=""))
}

write.csv(KKT.perc, file="KKT8.csv")
KKT.perc


load("b0_9.rda")
load("beta_9.rda")
load("lambda_9.rda")

dim(betas_9)
dim(betas_9)

kk = 2:5
KKT.perc = rep(NA, 4)
for(i in 1:length(kk))
{
  KKT.perc[i] = KKT(b0_9, betas_9, dat$y, dat$x, lambda_9, lambda2=l2, thr=10^(-kk[i]), 
                    qv=qvalue, loss = c("power"), print.out=F)
  write.csv(KKT.perc[i], file=paste("thr9_", kk[i], ".csv", sep=""))
}

write.csv(KKT.perc, file="KKT9.csv")
KKT.perc


load("b0_10.rda")
load("beta_10.rda")
load("lambda_10.rda")

dim(betas_10)
dim(betas_10)

kk = 2:5
KKT.perc = rep(NA, 4)
for(i in 1:length(kk))
{
  KKT.perc[i] = KKT(b0_10, betas_10, dat$y, dat$x, lambda_10, lambda2=l2, thr=10^(-kk[i]), 
                    qv=qvalue, loss = c("power"), print.out=F)
  write.csv(KKT.perc[i], file=paste("thr10_", kk[i], ".csv", sep=""))
}

write.csv(KKT.perc, file="KKT10.csv")

KKT.perc

stop1 = Sys.time()
difftime(stop1, start1, units="secs")
