##############################################################################
## 	Code for replicating Table 2: 
## 		Timings (in seconds) for WZZ and hubernet for n = 100, p = 5000 data.
##
## 	Note: we use WZZ (R code) by Wang et al. (2008) for the HHSVM,
## 		which is not public domain on CRAN.
## 		Hence you have to obtain this code to run it, 
##
## 	The source code "hhsvm.r" is available at
## 		http://www.stat.lsa.umich.edu/~jizhu/code/hhsvm/.
##
##  WZZ does not divide the loss by N, but with a factor delta for the loss. 
##    	So we must select the appropriate lambda1 and lambda2 values 
##  	for WZZ such that it produce the same solution as hubernet does.
##
source("source.R") # R code for generating the FHT data
source("hhsvm.r")  # load WZZ source code 
library(gcdnet) # load gcdnet package
# 	Note: need to make sure that 
#         gcdnet and hhsvm.r are installed on you system
require(MASS)   

nrep = 10 # the number of independent runs 
delta = 2 # delta for huberized hinge loss

n=100 # the number of observations
p=5000 # the number of predictors

##############################################################################
## lambda2 = 0 case
lambda2 = 0*n*delta
set.seed(44)
seeds=sample(1:1000,size=nrep*6)
tim=tim1=NULL
iii=0
for(ii in 1:nrep){
	for(rho in c(0,.1,.2,.5,.8,.95)){
		cat(rho,fill=T)
		iii=iii+1
		set.seed(seeds[iii])
		x=genx2(n,p,rho)
		y=genjerry(x,3)
		y=as.vector(1*(runif(n)< 1/(1+exp(-y))))
		y[y==0]=-1
		# run DrHSVM in WZZ
		tim=c(tim,unix.time(m1<-DrHSVM(x,y,lambda2,delta=delta,type="lasso",scale=F))[1])
		lambda = m1$lambda1/(delta*n)
		lambda = lambda[-length(lambda)]
		# run hubernet in gcdnet
		tim1=c(tim1,unix.time(m2<-gcdnet(method="hhsvm",x=x,y=y,eps=1e-8,lambda2=lambda2/(n*delta),delta=delta,standardize=F,lambda=lambda,maxit=3e7))[1])
}}
tim=matrix(tim,nrow=nrep,byrow=T)
tim1=matrix(tim1,nrow=nrep,byrow=T)

tim.m=apply(tim,2,mean)
tim1.m=apply(tim1,2,mean)

tims=round(rbind(tim.m,tim1.m),3)
dimnames(tims)=list(c("WZZ","hubernet"),as.character(c(0, .1,.2,.5,.8,.95)))
tims


##############################################################################
## lambda2 = 1e-4 case
lambda2 = 1e-4*n*delta
set.seed(44)
seeds=sample(1:1000,size=nrep*6)
tim=tim1=NULL
iii=0
for(ii in 1:nrep){
	for(rho in c(0,.1,.2,.5,.8,.95)){
		cat(rho,fill=T)
		iii=iii+1
		set.seed(seeds[iii])
		x=genx2(n,p,rho)
		y=genjerry(x,3)
		y=as.vector(1*(runif(n)< 1/(1+exp(-y))))
		y[y==0]=-1
		# run DrHSVM in WZZ
		tim=c(tim,unix.time(m1<-DrHSVM(x,y,lambda2,delta=delta,type="lasso",scale=F))[1])
		lambda = m1$lambda1/(delta*n)
		lambda = lambda[-length(lambda)]
		# run hubernet in gcdnet
		tim1=c(tim1,unix.time(m2<-gcdnet(method="hhsvm",x=x,y=y,eps=1e-8,lambda2=lambda2/(n*delta),delta=delta,standardize=F,lambda=lambda,maxit=3e7))[1])
}}
tim=matrix(tim,nrow=nrep,byrow=T)
tim1=matrix(tim1,nrow=nrep,byrow=T)

tim.m=apply(tim,2,mean)
tim1.m=apply(tim1,2,mean)

tims=round(rbind(tim.m,tim1.m),3)
dimnames(tims)=list(c("WZZ","hubernet"),as.character(c(0, .1,.2,.5,.8,.95)))
tims


##############################################################################
## lambda2 = 1e-2 case
lambda2 = 1e-2*n*delta
set.seed(44)
seeds=sample(1:1000,size=nrep*6)
tim=tim1=NULL
iii=0
for(ii in 1:nrep){
	for(rho in c(0,.1,.2,.5,.8,.95)){
		cat(rho,fill=T)
		iii=iii+1
		set.seed(seeds[iii])
		x=genx2(n,p,rho)
		y=genjerry(x,3)
		y=as.vector(1*(runif(n)< 1/(1+exp(-y))))
		y[y==0]=-1
		# run DrHSVM in WZZ
		tim=c(tim,unix.time(m1<-DrHSVM(x,y,lambda2,delta=delta,type="lasso",scale=F))[1])
		lambda = m1$lambda1/(delta*n)
		lambda = lambda[-length(lambda)]
		# run hubernet in gcdnet
		tim1=c(tim1,unix.time(m2<-gcdnet(method="hhsvm",x=x,y=y,eps=1e-8,lambda2=lambda2/(n*delta),delta=delta,standardize=F,lambda=lambda,maxit=3e7))[1])
}}
tim=matrix(tim,nrow=nrep,byrow=T)
tim1=matrix(tim1,nrow=nrep,byrow=T)

tim.m=apply(tim,2,mean)
tim1.m=apply(tim1,2,mean)

tims=round(rbind(tim.m,tim1.m),3)
dimnames(tims)=list(c("WZZ","hubernet"),as.character(c(0, .1,.2,.5,.8,.95)))
tims


##############################################################################
## lambda2 = 1 case
lambda2 = 1*n*delta
set.seed(44)
seeds=sample(1:1000,size=nrep*6)
tim=tim1=NULL
iii=0
for(ii in 1:nrep){
	for(rho in c(0,.1,.2,.5,.8,.95)){
		cat(rho,fill=T)
		iii=iii+1
		set.seed(seeds[iii])
		x=genx2(n,p,rho)
		y=genjerry(x,3)
		y=as.vector(1*(runif(n)< 1/(1+exp(-y))))
		y[y==0]=-1
		# run DrHSVM in WZZ
		tim=c(tim,unix.time(m1<-DrHSVM(x,y,lambda2,delta=delta,type="lasso",scale=F))[1])
		lambda = m1$lambda1/(delta*n)
		lambda = lambda[-length(lambda)]
		# run hubernet in gcdnet
		tim1=c(tim1,unix.time(m2<-gcdnet(method="hhsvm",x=x,y=y,eps=1e-8,lambda2=lambda2/(n*delta),delta=delta,standardize=F,lambda=lambda,maxit=3e7))[1])
}}
tim=matrix(tim,nrow=nrep,byrow=T)
tim1=matrix(tim1,nrow=nrep,byrow=T)

tim.m=apply(tim,2,mean)
tim1.m=apply(tim1,2,mean)

tims=round(rbind(tim.m,tim1.m),3)
dimnames(tims)=list(c("WZZ","hubernet"),as.character(c(0, .1,.2,.5,.8,.95)))
tims
