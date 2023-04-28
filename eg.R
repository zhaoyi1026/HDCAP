###########################################
# High-dimensional CAP
###########################################

library("mvtnorm")

rm(list=ls())

source("CAP_HD.R")

###########################################
# model parameters
load("par.RData")

n<-100
nT<-100

Tvec<-rep(nT,n) # number of time points for each subject
###########################################

#########################################
# generate data
set.seed(100)
X<-cbind(rep(1,n),rbinom(n,size=1,prob=0.5))

# Sigma matrix
Sigma<-array(NA,c(p,p,n))
delta<-matrix(NA,n,p)
i<-1
j<-1
for(i in 1:n)
{
  Delta<-matrix(0,p,p)
  for(j in 1:ncol(beta.mat))
  {
    if(beta.mat[2,j]==0)
    {
      Delta[j,j]<-exp(rnorm(1,mean=beta.mat[1,j],sd=beta.sd))
    }else
    {
      Delta[j,j]<-exp(t(X[i,])%*%beta.mat[,j])
    }
  }
  delta[i,]<-diag(Delta)
  Sigma[,,i]<-Gamma%*%Delta%*%t(Gamma)
}

# Y
set.seed(100)
Y<-list()
for(i in 1:n)
{
  Y[[i]]<-rmvnorm(n=Tvec[i],mean=rep(0,p),sigma=Sigma[,,i])
}
#########################################


###########################################
# method parameters
score.return<-TRUE
verbose<-TRUE

# bootstrap parameters
boot<-TRUE
sims<-500
boot.ci.type<-"perc"
conf.level<-0.95

# number of directions
DfD.thred<-2
###########################################

#########################################
# HD-CAP regression

# parameter estimation
re<-capReg(Y,X,stop.crt="DfD",DfD.thred=DfD.thred,cov.shrinkage=TRUE,score.return=score.return,verbose=verbose)

# to compare with the truth
# the 2nd and 4th dimension in Gamma are the components with nonzero beta (see beta.mat)
t(re$gamma)%*%Gamma
re$beta            # beta coefficient

# bootstrap inference for beta
re.boot<-vector("list",length=ncol(re$gamma))
for(jj in 1:ncol(re$gamma))
{
  re.boot[[jj]]<-cap_beta_boot(Y,X,gamma=re$gamma[,jj],cov.shrinkage=TRUE,boot=boot,sims=sims,boot.ci.type=boot.ci.type,conf.level=conf.level,verbose=TRUE)
}
for(jj in 1:ncol(re$gamma))
{
  print(re.boot[[jj]]$Inference)
}
#########################################
save.image("eg.RData")



