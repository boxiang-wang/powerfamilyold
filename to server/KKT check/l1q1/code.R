start1 = Sys.time()

l2 = 1
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
