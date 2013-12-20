rm(list=ls(all=TRUE))


#setwd("D:\\GitHub\\powerfamily")
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


# dyn.load("M_powerfamilyNET.dll")
dyn.load("M_powerfamilyNET.so")

args=(commandArgs(TRUE))
if(length(args)==0){
  print("No arguments supplied.")
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

#qv = 1
print(qv)
#FHT = FHTgen(n=100, p=5000, rho=0.8)
load("colon.rda")
x=colon.x
y=colon.y


nrep = 10
seed = sample((1:nrep) + 2001, nrep)
ans = matrix(0, 7, nrep)
time = matrix(0, 7, nrep)
for (j in 1:nrep){
  save(j, file=paste(j, "th.rda", sep=""))
  save(ans, file="anstemp.rda")
  save(time, file="timetemp.rda")
  i=1
  set.seed(seed[j])
  index = sample(1:nrow(x),as.integer(nrow(x)/3),replace=F) 
  test_x = x[index,]
  test_y = y[index]
  train_x = x[-index,]
  train_y = y[-index]
  for(lambda2 in c(1e-4,1e-3,1e-2,1e-1,1,5,10)){
    tim = system.time(cv<-cv.GCDpower(train_x, train_y, eps=1e-8, qv=qv,
                                      lambda2=lambda2, method="power",
                                      pred.loss="misclass", nfolds=5))[3]
    pre = predict(cv$GCDpower.fit, newx = test_x, s = cv$lambda.1se, type = "class" )
    error = (test_y != pre) 
    ans[i,j] = res = sum(error)/length(error)
    time[i,j] = tim 
    i = i + 1
  }}
(ans.avg = apply(ans,1,mean)*100)
(time.avg = apply(time,1,mean))

save(ans.avg, file="ans.rda")
save(time.avg, file="time.rda")

