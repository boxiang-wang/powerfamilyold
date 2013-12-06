rm(list=ls(all=TRUE))

#setwd("D:\\GitHub\\powerfamily")
require(Matrix)

# Main program
source("M_GCDpower.R")
source("O_utilities.R")
source("M_FHTgen.R")

dyn.load("M_powerfamilyintNET.so")
dyn.load("M_powerfamilyNET.so")


FHT = FHTgen(n=100, p=5000, rho=0.8)
dat = FHT
x = dat$x
y = dat$y

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

setwd("D:\\GitHub\\powerfamily\\to server\\timing\\t3")
load("avgtable.rda")
t3 = avg.time.table

setwd("D:\\GitHub\\powerfamily\\to server\\timing\\t5")
load("avgtable.rda")
t5 = avg.time.table

setwd("D:\\GitHub\\powerfamily\\to server\\timing\\t6")
load("avgtable.rda")
t6 = avg.time.table

t_n100p5Kr08 = (t1+t2)/2
library(xtable)
xtable(t_n100p5Kr08, digits=3)


t_n5Kp100r05 = (t3+t6)/2
library(xtable)
xtable(t_n5Kp100r05, digits=3)