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
source("M_cv.GCDpower.R")
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

print(qv)
#FHT = FHTgen(n=100, p=5000, rho=0.8)
load("data/colon.rda")
x = colon.x
y = colon.y
y = c(-1,1)[as.factor(y)]



nrep = 55
seed = sample((1:nrep) + 2001, nrep)
ans = matrix(0, 7, nrep)
time = matrix(0, 7, nrep)
nzo = matrix(0, 7, nrep)
for (j in 1:nrep){
  print(paste("q =", qv, "j =", j))
  write.csv(NA, file=paste("q_", qv, "_j_", j, ".csv", sep=""))
  i = 1
  set.seed(seed[j])
  index = sample(1:nrow(x),as.integer(nrow(x)/5),replace=F) 
  test_x = x[index,]
  test_y = y[index]
  train_x = x[-index,]
  train_y = y[-index]
  for(lambda2 in c(1e-4,1e-3,1e-2,1e-1,1,5,10)){
    tim = system.time(cv<-cv.GCDpower(train_x, train_y, eps=1e-8, qv=qv, delta=2,
                                      lambda2=lambda2, method="power",
                                      pred.loss="misclass", nfolds=5))[3]
    coef.cvs = coef(cv, s="lambda.1se")[-1,]
    nzov= length(coef.cvs[coef.cvs != 0])
    pre = predict(cv$GCDpower.fit, newx = test_x, s = cv$lambda.1se, type = "class" )
    error = (test_y != pre) 
    nzo[i,j] = nzov
    time[i,j] = tim
    ans[i,j] = res = mean(error)
    i = i + 1
  }
  save(ans, file=paste("q=", qv, "_ans.rda", sep=""))
  save(nzo, file=paste("q=", qv, "_nzo.rda", sep=""))
  save(time, file=paste("q=", qv, "_time.rda", sep=""))
}
(ans.avg = apply(ans,1,mean)*100)
(nzo.avg = apply(nzo,1,mean))
(time.avg = apply(time,1,mean))

save(ans.avg, file=paste("ans", qv, ".rda", sep=""))
save(nzo.avg, file=paste("nzo", qv, ".rda", sep=""))
save(time.avg, file=paste("time", qv, ".rda", sep=""))

