rm(list=ls(all=TRUE))

#setwd("D:\\GitHub\\powerfamily")
require(Matrix)

# Main program
source("M_GCDpower.R")
source("M_FHTgen.R")

#shell("Rcmd SHLIB M_powerfamilyNET.f90 M_powerfamilyintNET.f90 M_powerfamilyhalfNET.f90 O_auxiliary.f90 -o M_powerfamilyNET.so")

dyn.load("M_powerfamilyNET.dll")
dyn.load("M_powerfamilyNET.so")

set.seed(1234)
seeds = sample(1:1000, size=10*4*6)
FHT = FHTgen(n=100, p=5000, rho=0.8)
dat = FHT
x = dat$x
y = dat$y
seed.ind = 1
total.indep = 10

l2.list = c(0, 10^(-4), 10^(-2),1)
qv.list = c(0.25, 0.5, 1, 1.5, 2, 5)
avg.time.table = matrix(0, length(l2.list), length(qv.list))

for(indp in 1:total.indep)
{
  print(paste(indp, " th independent run.", sep=""))
  time.table = matrix(NA, length(l2.list), length(qv.list))
  for(i in 1:length(l2.list)) # row
  {
    for(j in 1:length(qv.list)) # column
    {
      
      l2 = l2.list[i]
      qv = qv.list[j]
      set.seed(seeds[seed.ind])
      FHT = FHTgen(n=100, p=5000, rho=0.8)
      dat = FHT
      x = dat$x
      y = dat$y
      seed.ind = seed.ind + 1
      start1 = Sys.time()
      m = gcdnetpower(x=x, y=y,
                      lambda2=l2, qv=qv, method="power",eps=1e-8, standardize=F)
      stop1 = Sys.time()
      time.table[i,j] = difftime(stop1, start1, units="secs")
    }
  }
  write.csv(time.table, file=paste("timetable_", indp, ".csv", sep=""))
  avg.time.table = avg.time.table + time.table
}


avg.time.table = avg.time.table / total.indep
write.csv(avg.time.table, file="avg.time.table.csv")
save(avg.time.table, file="avgtable.rda")




#start1 = Sys.time()
#m = gcdnetpower(x=dat$x, y=dat$y,
#                lambda2=1.5, qv=2, method="power",eps=1e-8, standardize=F)
#stop1 = Sys.time()
#difftime(stop1, start1, units="secs")

# combining outputs from server
setwd("D:\\GitHub\\powerfamily\\to server\\timing\\t1")
load("avgtable.rda")
t1 = avg.time.table

setwd("D:\\GitHub\\powerfamily\\to server\\timing\\t2")
load("avgtable.rda")
t2 = avg.time.table

setwd("D:\\GitHub\\powerfamily\\to server\\timing\\t11")
load("avgtable.rda")
t11 = avg.time.table

setwd("D:\\GitHub\\powerfamily\\to server\\timing\\t22")
load("avgtable.rda")
t22 = avg.time.table



t_n100p5Kr08 = (t1+t2)/2
library(xtable)
xtable(t_n100p5Kr08, digits=3)


setwd("D:\\GitHub\\powerfamily\\to server\\timing\\t3")
load("avgtable.rda")
t3 = avg.time.table

setwd("D:\\GitHub\\powerfamily\\to server\\timing\\t4")
load("avgtable.rda")
t4 = avg.time.table

setwd("D:\\GitHub\\powerfamily\\to server\\timing\\t33")
load("avgtable.rda")
t33 = avg.time.table

setwd("D:\\GitHub\\powerfamily\\to server\\timing\\t44")
load("avgtable.rda")
t44 = avg.time.table

t_n5Kp100r05 = (t3+t4)/2
library(xtable)
xtable(t_n5Kp100r05, digits=3)