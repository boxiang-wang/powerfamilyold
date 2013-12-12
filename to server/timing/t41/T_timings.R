rm(list=ls(all=TRUE))

#setwd("D:\\GitHub\\powerfamily")
require(Matrix)

# Main program
source("M_GCDpower.R")
source("O_utilities.R")
source("M_FHTgen.R")

dyn.load("M_powerfamilyNET.so")


FHT = FHTgen(n=5000, p=100, rho=0.5)
dat = FHT
x = dat$x
y = dat$y

total.indep = 5

l2.list = c(0, 10^(-4), 10^(-2),1)
qv.list = c(0.25, 0.5, 1, 2, 3, 5)
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


