#################################################
# High-dimensional Covariance regression
# linear shrinkage on the covariance matrix
#################################################

library("MASS")       # general inverse of a matrix
# library("nloptr")     # non-linear optimization
library("multigroup") # common PCA

library("glmnet")     # lasso package
#################################################
# soft-thresholding function
soft_thred<-function(mu,lambda)
{
  return(apply(cbind(abs(mu)-lambda,rep(0,length(mu))),1,function(x){return(max(x))})*sign(mu))
}

# standardized Frobenius norm
norm.F.std<-function(A1,A2=NULL)
{
  p<-nrow(A1)
    
  if(is.null(A2))
  {
    return(sqrt(sum(diag(A1%*%t(A1)))/p))
  }else
  {
    return(sum(diag(A1%*%t(A2)))/p)
  }
}
#################################################

#################################################
# objective function
obj.func<-function(Y,X,gamma,beta,lambda=0)
{
  n<-length(Y)
  p<-ncol(Y[[1]])
  Tvec<-rep(NA,n)
  
  q<-ncol(X)
  
  # Estimate covariance matrix for each subject
  Sigma<-array(NA,c(p,p,n))
  H<-matrix(0,p,p)
  for(i in 1:n)
  {
    Tvec[i]<-nrow(Y[[i]])
    
    Sigma[,,i]<-t(scale(Y[[i]],center=TRUE,scale=FALSE))%*%(scale(Y[[i]],center=TRUE,scale=FALSE))/nrow(Y[[i]])
    
    H<-H+Sigma[,,i]*Tvec[i]
  }
  H<-H/sum(Tvec)
  
  score<-apply(Sigma,3,function(x){return(t(gamma)%*%x%*%gamma)})
  
  # linear shrinkage on covariance matrix
  mu<-mean(exp(X%*%beta))/(t(gamma)%*%gamma)[1,1]
  delta2<-mean((score-mu*(t(gamma)%*%gamma)[1,1])^2)
  psi2<-min(mean((score-c(exp(X%*%beta)))^2),delta2)
  phi2<-delta2-psi2
  # phi2<-mean((mu*(t(gamma)%*%gamma)[1,1]-exp(X%*%beta))^2)
  rho1<-psi2*mu/delta2
  rho2<-phi2/delta2
  
  obj<-0
  for(i in 1:n)
  {
    Stmp<-rho1*diag(rep(1,p))+rho2*Sigma[,,i]
    obj<-obj+Tvec[i]*(t(X[i,])%*%beta+exp(-t(X[i,])%*%beta)*(t(gamma)%*%Stmp%*%gamma))[1,1]
  }
  
  pen<-2*lambda*sum(abs(gamma))
  
  re<-data.frame(obj=obj,pen=pen,sum=obj+pen)
  
  return(re)
}

# linear shrinkage of covariance matrix
cov.ls<-function(X)
{
  n<-nrow(X)
  p<-ncol(X)
  
  # demean of X
  X<-scale(X,center=TRUE,scale=FALSE)
  
  # sample covariance matrix
  S<-cov(X)*(n-1)/n
  
  Ip<-diag(rep(1,p))
  
  m<-norm.F.std(S,diag(rep(1,p)))
  d2<-(norm.F.std(S-m*diag(rep(1,p))))^2
  
  b2.bar<-mean(apply(X,1,function(x){return((norm.F.std(x%*%t(x)-S))^2)}))/n
  
  b2<-min(b2.bar,d2)
  
  a2<-d2-b2
  
  return(b2*m*Ip/d2+a2*S/d2)
}

# given gamma, estimate beta
cap_beta<-function(Y,X,gamma,cov.shrinkage=TRUE,max.itr=1000,tol=1e-4,trace=FALSE,score.return=TRUE)
{
  n<-length(Y)
  p<-ncol(Y[[1]])
  Tvec<-rep(NA,n)
  for(i in 1:n)
  {
    Tvec[i]<-nrow(Y[[i]])
  }
  
  q<-ncol(X)
  if(is.null(colnames(X))==TRUE)
  {
    colnames(X)<-paste0("X",1:q)
  }
  
  if(min(Tvec)-5<p)
  {
    cov.shrinkage<-TRUE
  }
  
  # Estimate covariance matrix for each subject
  Sigma<-array(NA,c(p,p,n))
  H<-matrix(0,p,p)
  for(i in 1:n)
  {
    if(cov.shrinkage)
    {
      Sigma[,,i]<-cov.ls(Y[[i]])
    }else
    {
      Sigma[,,i]<-cov(Y[[i]])*(Tvec[i]-1)/Tvec[i] 
    }
    
    H<-H+Sigma[,,i]*Tvec[i]
  }
  H<-H/sum(Tvec)
  
  score<-apply(Sigma,3,function(x){return(t(gamma)%*%x%*%gamma)})
  
  set.seed(100)
  beta0<-c(0,rnorm(q-1,mean=1,sd=1))
  # beta0<-rep(0,q)
  
  if(cov.shrinkage)
  {
    #-----------------------------------
    # estimate linear weights
    G<-diag(rep(1,p))
    mu<-mean(exp(X%*%beta0))/(t(gamma)%*%G%*%gamma)[1,1]
    delta2<-mean((score-mu*(t(gamma)%*%G%*%gamma)[1,1])^2)
    psi2<-min(mean((score-c(exp(X%*%beta0)))^2/Tvec),delta2)
    phi2<-delta2-psi2
    # psi2<-mean((score-c(exp(X%*%beta0)))^2/Tvec)
    # phi2<-mean((mu*(t(gamma)%*%G%*%gamma)[1,1]-exp(X%*%beta0))^2)
    rho1<-psi2*mu/(psi2+phi2)
    rho2<-phi2/(psi2+phi2)
    #-----------------------------------
  }

  if(trace)
  {
    beta.trace<-beta0
    
    if(cov.shrinkage)
    {
      rho.trace<-c(rho1,rho2) 
    }
    
    obj<-obj.func(Y=Y,X=X,gamma=gamma,beta=beta0,lambda=0)
  }
  
  s<-0
  diff<-100
  while(s<=max.itr&diff>tol)
  {
    s<-s+1
    
    # update beta
    #---------------------
    # likelihood method
    Q1<-matrix(0,q,q)
    Q2<-rep(0,q)
    for(i in 1:n)
    {
      if(cov.shrinkage)
      {
        Stmp<-rho1*G+rho2*Sigma[,,i] 
      }else
      {
        Stmp<-Sigma[,,i]
      }

      Q1<-Q1+(Tvec[i]*(t(gamma)%*%Stmp%*%gamma)[1,1]*exp(-t(X[i,])%*%beta0)[1,1])*(X[i,]%*%t(X[i,]))
      Q2<-Q2+Tvec[i]*(1-(t(gamma)%*%Stmp%*%gamma)[1,1]*(exp(-t(X[i,])%*%beta0)[1,1]))*X[i,]
    }
    # beta.new<-beta0-solve(Q1)%*%Q2
    beta.new<-beta0-ginv(Q1)%*%Q2
    #---------------------
    
    if(cov.shrinkage)
    {
      # estimate linear weights
      mu<-mean(exp(X%*%beta.new))/(t(gamma)%*%G%*%gamma)[1,1]
      delta2<-mean((score-mu*(t(gamma)%*%G%*%gamma)[1,1])^2)
      psi2<-min(mean((score-c(exp(X%*%beta.new)))^2/Tvec),delta2)
      phi2<-delta2-psi2
      # psi2<-mean((score-c(exp(X%*%beta.new)))^2/Tvec)
      # phi2<-mean((mu*(t(gamma)%*%G%*%gamma)[1,1]-exp(X%*%beta.new))^2)
      rho1<-psi2*mu/(psi2+phi2)
      rho2<-phi2/(psi2+phi2)
    }
    
    if(trace)
    {
      beta.trace<-cbind(beta.trace,beta.new)
      
      if(cov.shrinkage)
      {
        rho.trace<-cbind(rho.trace,c(rho1,rho2))
      }
      
      obj<-rbind(obj,obj.func(Y=Y,X=X,gamma=gamma,beta=beta.new,lambda=0))
    }
    
    diff<-max(abs(beta.new-beta0))
    
    beta0<-beta.new
    
    # print(c(diff,rho1,rho2))
  }
  
  if(cov.shrinkage)
  {
    score.sk<-rep(NA,n)
    for(i in 1:n)
    {
      Stmp<-rho1*G+rho2*Sigma[,,i] 
      score.sk[i]<-t(gamma)%*%Stmp%*%gamma
    }
  }
  
  if(trace)
  {
    rownames(beta.trace)<-colnames(X)
    if(cov.shrinkage)
    {
      rownames(rho.trace)<-c("rho1","rho2")
    }
    if(score.return)
    {
      if(cov.shrinkage)
      {
        re<-list(gamma=c(gamma),beta=c(beta.new),convergence=(s<max.itr),shrinkage=data.frame(rho1=rho1,rho2=rho2,mu=mu,phi2=phi2,psi2=psi2,delta2=delta2),score=score,score.shrinkage=score.sk,
                 beta.trace=beta.trace,rho.trace=rho.trace,obj=obj)  
      }else
      {
        re<-list(gamma=c(gamma),beta=c(beta.new),convergence=(s<max.itr),score=score,beta.trace=beta.trace,obj=obj)
      }

    }else
    {
      if(cov.shrinkage)
      {
        re<-list(gamma=c(gamma),beta=c(beta.new),convergence=(s<max.itr),shrinkage=data.frame(rho1=rho1,rho2=rho2,mu=mu,phi2=phi2,psi2=psi2,delta2=delta2),beta.trace=beta.trace,rho.trace=rho.trace,obj=obj)
      }else
      {
        re<-list(gamma=c(gamma),beta=c(beta.new),convergence=(s<max.itr),beta.trace=beta.trace,obj=obj)
      }
    }
  }else
  {
    if(score.return)
    {
      if(cov.shrinkage)
      {
        re<-list(gamma=c(gamma),beta=c(beta.new),convergence=(s<max.itr),shrinkage=data.frame(rho1=rho1,rho2=rho2,mu=mu,phi2=phi2,psi2=psi2,delta2=delta2),score=score,score.shrinkage=score.sk)  
      }else
      {
        re<-list(gamma=c(gamma),beta=c(beta.new),convergence=(s<max.itr),score=score) 
      }
    }else
    {
      if(cov.shrinkage)
      {
        re<-list(gamma=c(gamma),beta=c(beta.new),convergence=(s<max.itr),shrinkage=data.frame(rho1=rho1,rho2=rho2,mu=mu,phi2=phi2,psi2=psi2,delta2=delta2)) 
      }else
      {
        re<-list(gamma=c(gamma),beta=c(beta.new),convergence=(s<max.itr))
      }
    }
  }
  
  return(re)
}
#################################################

#################################################
# estimate both gamma and beta
# finding first direction
gamma.solve<-function(A,H)
{
  p<-ncol(H)
  
  H.svd<-svd(H)
  H.d.sqrt<-diag(sqrt(H.svd$d))
  H.d.sqrt.inv<-diag(1/sqrt(H.svd$d))
  H.sqrt.inv<-H.svd$u%*%H.d.sqrt.inv%*%t(H.svd$u)
  
  #---------------------------------------------------
  # svd decomposition method
  eigen.tmp<-eigen(H.d.sqrt.inv%*%t(H.svd$u)%*%A%*%H.svd$u%*%H.d.sqrt.inv)
  eigen.tmp.vec<-Re(eigen.tmp$vectors)
  
  # obj<-rep(NA,ncol(eigen.tmp$vectors))
  # for(j in 1:ncol(eigen.tmp$vectors))
  # {
  #   otmp<-H.svd$u%*%H.d.sqrt.inv%*%eigen.tmp.vec[,j]
  #   obj[j]<-t(otmp)%*%A%*%otmp
  # }
  # re<-H.svd$u%*%H.d.sqrt.inv%*%eigen.tmp.vec[,which.min(obj)]
  re<-H.svd$u%*%H.d.sqrt.inv%*%eigen.tmp.vec[,p]
  #---------------------------------------------------
  
  #---------------------------------------------------
  # eigenvector of A with respect to H
  # eigen.tmp<-eigen(H.sqrt.inv%*%A%*%H.sqrt.inv)
  # eigen.tmp.vec<-Re(eigen.tmp$vectors)
  # 
  # obj<-rep(NA,ncol(eigen.tmp$vectors))
  # for(j in 1:ncol(eigen.tmp$vectors))
  # {
  #   otmp<-H.sqrt.inv%*%eigen.tmp.vec[,j]
  #   obj[j]<-t(otmp)%*%A%*%otmp
  # }
  # re<-H.sqrt.inv%*%eigen.tmp.vec[,p]
  #---------------------------------------------------
  
  return(re)
}

# gamma thresholding function
gamma.thred<-function(gamma,lambda,Y=NULL,method=c("HT","Reg"))
{
  # standardize gamma first
  gamma<-gamma/sqrt(sum(gamma^2))
  
  # hard-thresholding of lambda
  if(method=="HT")
  {
    gamma.th<-gamma
    gamma.th[abs(gamma.th)<lambda]<-0
    
    if(sqrt(sum(gamma.th^2))!=0)
    {
      gamma.th<-gamma.th/sqrt(sum(gamma.th^2))
    }else
    {
      while(sqrt(sum(gamma.th^2))==0)
      {
        lambda<-lambda/2
        
        gamma.th<-gamma
        gamma.th[abs(gamma.th)<lambda]<-0
      }
      
      gamma.th<-gamma.th/sqrt(sum(gamma.th^2))
    }
  }
  
  # regression based method
  if(method=="Reg")
  {
    # concatenate all observations
    Ytmp<-NULL
    Ztmp<-NULL
    for(i in 1:length(Y))
    {
      Ytmp<-rbind(Ytmp,Y[[i]])
      Ztmp<-c(Ztmp,Y[[i]]%*%gamma)
    }
    out.glmnet<-glmnet(x=Ytmp,y=Ztmp,family="gaussian",lambda=lambda)
    
    gamma.th<-as.vector(out.glmnet$beta)
    
    if(sqrt(sum(gamma.th^2))!=0)
    {
      gamma.th<-gamma.th/sqrt(sum(gamma.th^2))
    }else
    {
      while(sqrt(sum(gamma.th^2))==0)
      {
        lambda<-lambda/2
        
        out.glmnet<-glmnet(x=Ytmp,y=Ztmp,family="gaussian",lambda=lambda)
        gamma.th<-as.vector(out.glmnet$beta)
      }
      
      gamma.th<-gamma.th/sqrt(sum(gamma.th^2))
    }
  }
  
  re<-list(method=method,gamma=gamma,gamma.thred=gamma.th,lambda=lambda)
}

# estimate both gamma and beta
cap_D1<-function(Y,X,method=c("CAP"),lambda=0,method.thred=c("HT","Reg"),cov.shrinkage=FALSE,max.itr=1000,tol=1e-4,score.return=TRUE,trace=FALSE,gamma0=NULL)
{
  n<-length(Y)
  p<-ncol(Y[[1]])
  Tvec<-rep(NA,n)
  for(i in 1:n)
  {
    Tvec[i]<-nrow(Y[[i]])
  }
  
  q<-ncol(X)
  if(is.null(colnames(X))==TRUE)
  {
    colnames(X)<-paste0("X",1:q)
  }
  
  if(min(Tvec)-5<p)
  {
    cov.shrinkage<-TRUE
  }
  
  # Estimate covariance matrix for each subject
  Sigma<-array(NA,c(p,p,n))
  H<-matrix(0,p,p)
  for(i in 1:n)
  {
    if(cov.shrinkage)
    {
      Sigma[,,i]<-cov.ls(Y[[i]])
    }else
    {
      Sigma[,,i]<-cov(Y[[i]])*(Tvec[i]-1)/Tvec[i] 
    }
    
    H<-H+Sigma[,,i]*Tvec[i]
  }
  H<-H/sum(Tvec)
  
  if(method[1]=="CAP")
  {
    set.seed(100)
    beta0<-c(0,rnorm(q-1,mean=1,sd=1))
    # beta0<-rep(0,q)
    if(is.null(gamma0))
    {
      gamma0<-rep(1/sqrt(p),p) 
    }
    
    score<-apply(Sigma,3,function(x){return(t(gamma0)%*%x%*%gamma0)})
    
    if(cov.shrinkage)
    {
      #---------------------------------
      # estimate linear weights
      G<-diag(rep(1,p))
      mu<-mean(exp(X%*%beta0))/(t(gamma0)%*%G%*%gamma0)[1,1]
      delta2<-mean((score-mu*(t(gamma0)%*%G%*%gamma0)[1,1])^2)
      psi2<-min(mean((score-c(exp(X%*%beta0)))^2/Tvec),delta2)
      phi2<-delta2-psi2
      # psi2<-mean((score-c(exp(X%*%beta0)))^2/Tvec)
      # phi2<-mean((mu*(t(gamma0)%*%G%*%gamma0)[1,1]-exp(X%*%beta0))^2)
      rho1<-psi2*mu/(psi2+phi2)
      rho2<-phi2/(psi2+phi2)
      #---------------------------------
    }
    
    if(trace)
    {
      beta.trace<-NULL
      gamma.trace<-NULL
      rho.trace<-NULL
      
      obj<-NULL
    }
    
    s<-0
    diff<-100
    while(s<=max.itr&diff>tol)
    {
      s<-s+1
      
      #--------------------------------
      # update beta
      Q1<-matrix(0,q,q)
      Q2<-rep(0,q)
      for(i in 1:n)
      {
        if(cov.shrinkage)
        {
          Stmp<-rho1*G+rho2*Sigma[,,i]  
        }else
        {
          Stmp<-Sigma[,,i]
        }
        
        Q1<-Q1+(Tvec[i]*(t(gamma0)%*%Stmp%*%gamma0)[1,1]*exp(-t(X[i,])%*%beta0)[1,1])*(X[i,]%*%t(X[i,]))
        Q2<-Q2+Tvec[i]*(1-(t(gamma0)%*%Stmp%*%gamma0)[1,1]*(exp(-t(X[i,])%*%beta0)[1,1]))*X[i,]
      }
      beta.new<-beta0-ginv(Q1)%*%Q2
      #--------------------------------
      
      if(cov.shrinkage)
      {
        #--------------------------------
        # estimate linear weights
        mu<-mean(exp(X%*%beta.new))/(t(gamma0)%*%G%*%gamma0)[1,1]
        delta2<-mean((score-mu*(t(gamma0)%*%G%*%gamma0)[1,1])^2)
        psi2<-min(mean((score-c(exp(X%*%beta.new)))^2/Tvec),delta2)
        phi2<-delta2-psi2
        # psi2<-mean((score-c(exp(X%*%beta.new)))^2/Tvec)
        # phi2<-mean((mu*(t(gamma0)%*%G%*%gamma0)[1,1]-exp(X%*%beta.new))^2)
        rho1<-psi2*mu/(psi2+phi2)
        rho2<-phi2/(psi2+phi2)
        #--------------------------------
      }
      
      #--------------------------------
      # update gamma
      S1<-matrix(0,p,p)
      for(i in 1:n)
      {
        if(cov.shrinkage)
        {
          Stmp<-rho1*G+rho2*Sigma[,,i]                                             # shrinkage covariance matrix
        }else
        {
          Stmp<-Sigma[,,i]
        }
        
        S1<-S1+Tvec[i]*exp(-(t(X[i,])%*%beta.new)[1,1])*Stmp
      }
      # gamma.new<-ginv(S1+nu2*H+diag(rep(rho,p)))%*%(nu1+rho*alpha0)            # not a good solution
      # gamma.new<-gamma.solve.alpha(S1,H,alpha0,nu1,rho)                        # solution also taken into consideration of alpha
      gamma.new<-gamma.solve(S1,H)
      # gamma.new<-gamma.new/sqrt(sum(gamma.new^2))
      # t(gamma.new)%*%H%*%gamma.new
      #--------------------------------
      
      if(trace)
      {
        beta.trace<-cbind(beta.trace,beta.new)
        gamma.trace<-cbind(gamma.trace,gamma.new)
        
        if(cov.shrinkage)
        {
          rho.trace<-cbind(rho.trace,c(rho1,rho2)) 
        }
        
        obj<-rbind(obj,obj.func(Y=Y,X=X,gamma=gamma.new,beta=beta.new,lambda=0))
      }
      
      diff<-max(c(max(abs(gamma.new-gamma0)),max(abs(beta.new-beta0))))
      
      gamma0<-gamma.new
      beta0<-beta.new
      
      # print(c(diff,rho1,rho2))
    }
    
    gamma.new<-c(gamma.new)/sqrt(sum(gamma.new^2))
    if(gamma.new[which.max(abs(gamma.new))]<0)
    {
      gamma.new<--gamma.new
    }
    out.thred<-gamma.thred(gamma.new,lambda,Y=Y,method=method.thred)
    
    gamma.est<-out.thred$gamma.thred
    
    beta.est<-cap_beta(Y,X,gamma.est,cov.shrinkage=cov.shrinkage,max.itr=max.itr,tol=tol,trace=FALSE,score.return=FALSE)$beta
    
    score<-apply(Sigma,3,function(x){return(t(gamma.est)%*%x%*%gamma.est)})
    if(cov.shrinkage)
    {
      score.sk<-apply(Sigma,3,function(x){return(t(gamma.est)%*%(rho1*G+rho2*x)%*%gamma.est)}) 
    }
    
    if(trace)
    {
      rownames(beta.trace)<-colnames(X)
      if(cov.shrinkage)
      {
        rownames(rho.trace)<-c("rho1","rho2") 
      }
      if(score.return)
      {
        if(cov.shrinkage)
        {
          re<-list(lambda=out.thred$lambda,gamma=gamma.est,beta=beta.est,convergence=(s<max.itr),shrinkage=data.frame(rho1=rho1,rho2=rho2,mu=mu,phi2=phi2,psi2=psi2,delta2=delta2),
                   score=score,score.shrinkage=score.sk,gamma.unsparse=gamma.new,gamma.trace=gamma.trace,beta.trace=beta.trace,rho.trace=rho.trace,obj=obj)
        }else
        {
          re<-list(lambda=out.thred$lambda,gamma=gamma.est,beta=beta.est,convergence=(s<max.itr),score=score,gamma.unsparse=gamma.new,gamma.trace=gamma.trace,beta.trace=beta.trace,rho.trace=rho.trace,obj=obj)
        }
      }else
      {
        if(cov.shrinkage)
        {
          re<-list(lambda=out.thred$lambda,gamma=gamma.est,beta=beta.est,convergence=(s<max.itr),shrinkage=data.frame(rho1=rho1,rho2=rho2,mu=mu,phi2=phi2,psi2=psi2,delta2=delta2),
                   gamma.unsparse=gamma.new,gamma.trace=gamma.trace,beta.trace=beta.trace,rho.trace=rho.trace,obj=obj)
        }else
        {
          re<-list(lambda=out.thred$lambda,gamma=gamma.est,beta=beta.est,convergence=(s<max.itr),gamma.unsparse=gamma.new,gamma.trace=gamma.trace,beta.trace=beta.trace,rho.trace=rho.trace,obj=obj)
        }
      }
      
    }else
    {
      if(score.return)
      {
        if(cov.shrinkage)
        {
          re<-list(lambda=out.thred$lambda,gamma=gamma.est,beta=beta.est,convergence=(s<max.itr),shrinkage=data.frame(rho1=rho1,rho2=rho2,mu=mu,phi2=phi2,psi2=psi2,delta2=delta2),
                   score=score,score.shrinkage=score.sk,gamma.unsparse=gamma.new)
        }else
        {
          re<-list(lambda=out.thred$lambda,gamma=gamma.est,beta=beta.est,convergence=(s<max.itr),score=score,gamma.unsparse=gamma.new)
        }
      }else
      {
        if(cov.shrinkage)
        {
          re<-list(lambda=out.thred$lambda,gamma=gamma.est,beta=beta.est,convergence=(s<max.itr),shrinkage=data.frame(rho1=rho1,rho2=rho2,mu=mu,phi2=phi2,psi2=psi2,delta2=delta2),gamma.unsparse=gamma.new)
        }else
        {
          re<-list(lambda=out.thred$lambda,gamma=gamma.est,beta=beta.est,convergence=(s<max.itr),gamma.unsparse=gamma.new)
        }
      }
    }
  }
  
  return(re)
}

# try several initial value of gamma and optimize over the objective function
cap_D1_opt<-function(Y,X,method=c("CAP"),lambda=0,method.thred=c("HT","Reg"),cov.shrinkage=FALSE,max.itr=1000,tol=1e-4,score.return=TRUE,trace=FALSE,gamma0.mat=NULL,ninitial=NULL,seed=500)
{
  n<-length(Y)
  p<-ncol(Y[[1]])
  Tvec<-rep(NA,n)
  for(i in 1:n)
  {
    Tvec[i]<-nrow(Y[[i]])
  }
  
  q<-ncol(X)
  if(is.null(colnames(X))==TRUE)
  {
    colnames(X)<-paste0("X",1:q)
  }
  
  if(min(Tvec)-5<p)
  {
    cov.shrinkage<-TRUE
  }
  
  # Estimate covariance matrix for each subject
  Sigma<-array(NA,c(p,p,n))
  H<-matrix(0,p,p)
  for(i in 1:n)
  {
    if(cov.shrinkage)
    {
      Sigma[,,i]<-cov.ls(Y[[i]])
    }else
    {
      Sigma[,,i]<-cov(Y[[i]])*(Tvec[i]-1)/Tvec[i] 
    }
    
    H<-H+Sigma[,,i]*Tvec[i]
  }
  H<-H/sum(Tvec)
  
  if(method[1]=="CAP")
  {
    # set initial values
    if(is.null(gamma0.mat))
    {
      #--------------------------------
      # gamma0.mat<-matrix(NA,p,p+1+5)
      # for(j in 1:p)
      # {
      #   gamma0.mat[,j]<-rep(0,p)
      #   gamma0.mat[j,j]<-1
      # }
      # gamma0.mat[,p+1]<-rep(1,p)/sqrt(sum(rep(1,p)^2))
      # 
      # set.seed(500)
      # gamma.tmp<-matrix(rnorm(5*p,mean=0,sd=1),nrow=p)
      # gamma0.mat[,(p+2):(p+1+5)]<-apply(gamma.tmp,2,function(x){return(x/sqrt(sum(x^2)))})
      #--------------------------------
      
      #--------------------------------
      # H eigenvalue
      # H.eigen<-eigen(H)
      # if(p>20)
      # {
      #   # set.seed(500)
      #   # gamma.tmp<-matrix(rnorm(25*p,mean=0,sd=1),nrow=p)
      #   # gamma0.mat<-apply(gamma.tmp,2,function(x){return(x/sqrt(sum(x^2)))})
      #   
      #   gamma0.mat<-H.eigen$vectors[,1:max(5,(sum(cumsum(H.eigen$values)/sum(H.eigen$values)<0.85)+1))]
      # }else
      # {
      #   gamma0.mat<-H.eigen$vectors
      # }
      
      set.seed(seed)
      gamma.tmp<-matrix(rnorm((p+1+5)*p,mean=0,sd=1),nrow=p)
      gamma0.mat<-apply(gamma.tmp,2,function(x){return(x/sqrt(sum(x^2)))})
      #--------------------------------
    }
    if(is.null(ninitial))
    {
      ninitial<-min(ncol(gamma0.mat),10)
    }else
    {
      if(ninitial>ncol(gamma0.mat))
      {
        ninitial<-ncol(gamma0.mat)
      }
    }
    set.seed(seed)
    gamma0.mat<-matrix(gamma0.mat[,sort(sample(1:ncol(gamma0.mat),ninitial,replace=FALSE))],ncol=ninitial)
    # gamma0.mat<-matrix(gamma0.mat[,1:ninitial],ncol=ninitial)
    
    re.tmp<-vector("list",ncol(gamma0.mat))
    obj<-rep(NA,ncol(gamma0.mat))
    for(kk in 1:ncol(gamma0.mat))
    {
      try(re.tmp[[kk]]<-cap_D1(Y,X,method=method[1],lambda=lambda,method.thred=method.thred[1],cov.shrinkage=cov.shrinkage,
                               max.itr=max.itr,tol=tol,score.return=score.return,trace=trace,gamma0=gamma0.mat[,kk]))
      
      if(is.null(re.tmp[[kk]])==FALSE)
      {
        gamma.unscale<-re.tmp[[kk]]$gamma/sqrt(t(re.tmp[[kk]]$gamma)%*%H%*%re.tmp[[kk]]$gamma)[1,1]
        
        try(beta.tmp<-cap_beta(Y,X,gamma.unscale,cov.shrinkage=cov.shrinkage,max.itr=max.itr,tol=tol,trace=FALSE,score.return=FALSE))
        try(obj[kk]<-obj.func(Y,X,gamma.unscale,beta.tmp$beta,lambda=0)$obj)
      }
    }
    opt.idx<-which.min(obj)
    re<-re.tmp[[opt.idx]]
  }
  
  return(re)
}
#################################################

#################################################
# second and higher direction
# CAP.OC: orthogonal constraint
cap_Dk<-function(Y,X,Phi0=NULL,method=c("CAP"),OC=FALSE,lambda=0,method.thred=c("HT","Reg"),cov.shrinkage=FALSE,max.itr=1000,tol=1e-4,score.return=TRUE,trace=FALSE,gamma0.mat=NULL,ninitial=NULL,seed=500)
{
  if(is.null(Phi0))
  {
    return(cap_D1_opt(Y,X,method=method[1],cov.shrinkage=cov.shrinkage,max.itr=max.itr,tol=tol,trace=trace,score.return=score.return,gamma0.mat=gamma0.mat,ninitial=ninitial,seed=seed))
  }else
  {
    n<-length(Y)
    q<-ncol(X)
    p<-ncol(Y[[1]])
    
    Tvec<-sapply(Y,nrow)
    
    if(min(Tvec)-5<p)
    {
      cov.shrinkage<-TRUE
    }
    
    p0<-ncol(Phi0)
    # estimate beta0
    beta0<-rep(NA,p0)
    for(j in 1:p0)
    {
      beta0[j]<-cap_beta(Y,X,gamma=Phi0[,j],cov.shrinkage=cov.shrinkage,max.itr=max.itr,tol=tol,trace=FALSE)$beta[1]
    }
    Ytmp<-vector("list",length=n)
    for(i in 1:n)
    {
      Y2tmp<-Y[[i]]-Y[[i]]%*%(Phi0%*%t(Phi0))
      ntmp<-nrow(Y2tmp)
      
      if(cov.shrinkage)
      {
        #----------------------------------
        # Y.svd.tmp<-svd(Y[[i]])
        # Y2tmp.svd<-svd(Y2tmp)
        # 
        # Y2tmp.eigen<-eigen(cov(Y2tmp))
        # ev.tmp<-Y2tmp.eigen$values
        # 
        # ev.new<-c(sqrt(ev.tmp[1:(ntmp-p0)]*ntmp),sqrt(exp(beta0)*ntmp))
        # 
        # Ytmp[[i]]<-Y2tmp.svd$u%*%diag(ev.new)%*%t(Y2tmp.eigen$vectors[,1:ntmp])
        #----------------------------------
        
        #----------------------------------
        # no need to add back the intercept using shrinkage method
        Ytmp[[i]]<-Y2tmp
        #----------------------------------
      }else
      {
        #----------------------------------
        Y2tmp.svd<-svd(Y2tmp)
        Ytmp[[i]]<-Y2tmp.svd$u%*%diag(c(Y2tmp.svd$d[1:(p-p0)],sqrt(exp(beta0)*ntmp)))%*%t(Y2tmp.svd$v)
        #----------------------------------
        
        #----------------------------------
        # Ytmp[[i]]<-Y2tmp
        #----------------------------------
      }
    }
    
    if(method[1]=="CAP")
    {
      if(OC==FALSE)
      {
        re.tmp<-cap_D1_opt(Ytmp,X,method=method[1],lambda=lambda,method.thred=method.thred[1],cov.shrinkage=cov.shrinkage,
                           max.itr=max.itr,tol=tol,score.return=score.return,trace=FALSE,gamma0.mat=gamma0.mat,ninitial=ninitial,seed=seed)
      }
    }else
    {
      re.tmp<-cap_D1_opt(Ytmp,X,method=method[1],lambda=lambda,method.thred=method.thred[1],cov.shrinkage=cov.shrinkage,
                         max.itr=max.itr,tol=tol,score.return=score.return,trace=FALSE,gamma0.mat=gamma0.mat,ninitial=ninitial,seed=seed) 
    }
    
    re<-re.tmp
    
    re$orthogonal<-data.frame(gamma=c(t(re.tmp$gamma)%*%Phi0),gamma.unsparse=c(t(re.tmp$gamma.unsparse)%*%Phi0))
    
    return(re)
  }
}
#################################################

#################################################
# level of diagonalization
diag.level<-function(Y,Phi)
{
  if(is.null(ncol(Phi))|ncol(Phi)==1)
  {
    stop("dimension of Phi is less than 2")
  }else
  {
    n<-length(Y)
    p<-ncol(Y[[1]])
    
    Tvec<-rep(NA,n)
    
    ps<-ncol(Phi)
    
    dl.sub<-matrix(NA,n,ps)
    colnames(dl.sub)<-paste0("Dim",1:ps)
    dl.sub[,1]<-1
    for(i in 1:n)
    {
      cov.tmp<-cov(Y[[i]])
      Tvec[i]<-nrow(Y[[i]])
      
      for(j in 2:ps)
      {
        phi.tmp<-Phi[,1:j]
        mat.tmp<-t(phi.tmp)%*%cov.tmp%*%phi.tmp
        dl.sub[i,j]<-det(diag(diag(mat.tmp)))/det(mat.tmp)
      }
    }
    
    pmean<-apply(dl.sub,2,function(y){return(prod(apply(cbind(y,Tvec),1,function(x){return(x[1]^(x[2]/sum(Tvec)))})))})
    
    re<-list(avg.level=pmean,sub.level=dl.sub)
    return(re)
  }
}
#################################################

#################################################
# finding first k directions
capReg<-function(Y,X,stop.crt=c("nD","DfD"),nD=NULL,DfD.thred=5,method=c("CAP"),OC=FALSE,lambda=0,method.thred=c("HT","Reg"),cov.shrinkage=FALSE,max.itr=1000,tol=1e-4,score.return=TRUE,trace=FALSE,
                 gamma0.mat=NULL,ninitial=NULL,seed=500,verbose=TRUE)
{
  # Y: outcome list
  # X: covariates
  # stop.crt: stopping criterion, nD=# of directions, DfD=DfD threshold
  
  if(stop.crt[1]=="nD"&is.null(nD))
  {
    stop.crt<-"DfD"
  }
  
  n<-length(Y)
  q<-ncol(X)
  p<-ncol(Y[[1]])
  
  if(is.null(colnames(X)))
  {
    colnames(X)<-paste0("X",1:q)
  }
  
  if(method[1]=="CAP")
  {
    #--------------------------------------------
    # First direction
    tm1<-system.time(re1<-cap_D1_opt(Y,X,method=method[1],lambda=lambda,method.thred=method.thred,cov.shrinkage=cov.shrinkage,
                                     max.itr=max.itr,tol=tol,score.return=score.return,trace=FALSE,gamma0.mat=gamma0.mat,ninitial=ninitial,seed=seed))
    
    Phi.est<-matrix(re1$gamma,ncol=1)
    beta.est<-matrix(re1$beta,ncol=1)
    Phi.est.unsparse<-matrix(re1$gamma.unsparse,ncol=1)
    
    cp.time<-matrix(as.numeric(tm1[1:3]),ncol=1)
    rownames(cp.time)<-c("user","system","elapsed")
    
    if(cov.shrinkage)
    {
      sk.out<-matrix(re1$shrinkage,nrow=1)
      colnames(sk.out)<-colnames(re1$shrinkage)
    }
    
    lambda.out<-re1$lambda
    
    if(verbose)
    {
      print(paste0("Component ",ncol(Phi.est)))
    }
    #--------------------------------------------
    
    if(stop.crt[1]=="nD")
    {
      if(score.return)
      {
        score<-matrix(re1$score,ncol=1)
        if(cov.shrinkage)
        {
          score.sk<-matrix(re1$score.shrinkage,ncol=1) 
        }
      }
      
      if(nD>1)
      {
        Phi.est<-matrix(Phi.est,ncol=1)
        for(j in 2:nD)
        {
          re.tmp<-NULL
          try(tm.tmp<-system.time(re.tmp<-cap_Dk(Y,X,Phi0=Phi.est,method=method,OC=OC,lambda=lambda,method.thred=method.thred,cov.shrinkage=cov.shrinkage,
                                                 max.itr=max.itr,tol=tol,score.return=score.return,trace=FALSE,gamma0.mat=gamma0.mat,ninitial=ninitial,seed=seed)))
          
          if(is.null(re.tmp)==FALSE)
          {
            Phi.est<-cbind(Phi.est,re.tmp$gamma)
            beta.est<-cbind(beta.est,re.tmp$beta)
            Phi.est.unsparse<-cbind(Phi.est.unsparse,re.tmp$gamma.unsparse)
            
            cp.time<-cbind(cp.time,as.numeric(tm.tmp[1:3]))
            
            if(cov.shrinkage)
            {
              sk.out<-rbind(sk.out,matrix(re.tmp$shrinkage,nrow=1))
            }
            
            lambda.out<-c(lambda.out,re.tmp$lambda)
            
            if(verbose)
            {
              print(paste0("Component ",ncol(Phi.est)))
            }
            
            if(score.return)
            {
              score<-cbind(score,re.tmp$score)  
              if(cov.shrinkage)
              {
                score.sk<-cbind(score.sk,re.tmp$score.shrinkage) 
              }
            }
          }else
          {
            break
          }
        }
      }
      
      colnames(Phi.est)=colnames(beta.est)=colnames(Phi.est.unsparse)<-paste0("D",1:ncol(Phi.est))
      if(cov.shrinkage)
      {
        rownames(sk.out)<-paste0("D",1:ncol(Phi.est)) 
      }
      rownames(Phi.est)=rownames(Phi.est.unsparse)<-paste0("V",1:p)
      rownames(beta.est)<-colnames(X)
      
      cp.time<-cbind(cp.time,apply(cp.time,1,sum))
      colnames(cp.time)<-c(paste0("D",1:ncol(Phi.est)),"Total")
      
      if(ncol(Phi.est)>1)
      {
        DfD<-diag.level(Y,Phi.est)
        DfD.unsparse<-diag.level(Y,Phi.est.unsparse)
      }else
      {
        DfD<-list(avg.level=1,sub.level=rep(1,n))
        DfD.unsparse<-list(avg.level=1,sub.level=rep(1,n))
      }
    }
    if(stop.crt[1]=="DfD")
    {
      nD<-1
      if(score.return)
      {
        score<-matrix(re1$score,ncol=1)
        if(cov.shrinkage)
        {
          score.sk<-matrix(re1$score.shrinkage,ncol=1) 
        }
      }
      
      DfD.tmp<-1
      while(DfD.tmp<DfD.thred)
      {
        re.tmp<-NULL
        try(tm.tmp<-system.time(re.tmp<-cap_Dk(Y,X,Phi0=Phi.est,method=method,OC=OC,lambda=lambda,method.thred=method.thred,cov.shrinkage=cov.shrinkage,
                                               max.itr=max.itr,tol=tol,score.return=score.return,trace=FALSE,gamma0.mat=gamma0.mat,ninitial=ninitial,seed=seed)))
        
        if(is.null(re.tmp)==FALSE)
        {
          nD<-nD+1
          
          DfD<-diag.level(Y,cbind(Phi.est,re.tmp$gamma))
          DfD.tmp<-DfD$avg.level[nD]
          
          if(DfD.tmp<DfD.thred)
          {
            Phi.est<-cbind(Phi.est,re.tmp$gamma)
            beta.est<-cbind(beta.est,re.tmp$beta)
            Phi.est.unsparse<-cbind(Phi.est.unsparse,re.tmp$gamma.unsparse)
            
            cp.time<-cbind(cp.time,as.numeric(tm.tmp[1:3]))
            
            if(cov.shrinkage)
            {
              sk.out<-rbind(sk.out,matrix(re.tmp$shrinkage,nrow=1))
            }
            
            lambda.out<-c(lambda.out,re.tmp$lambda)
            
            if(verbose)
            {
              print(paste0("Component ",ncol(Phi.est)))
            }
            
            if(score.return)
            {
              score<-cbind(score,re.tmp$score)  
              if(cov.shrinkage)
              {
                score.sk<-cbind(score.sk,re.tmp$score.shrinkage) 
              }
            }
          }
        }else
        {
          break
        }
      }
      
      colnames(Phi.est)=colnames(beta.est)=colnames(Phi.est.unsparse)<-paste0("D",1:ncol(Phi.est))
      if(cov.shrinkage)
      {
        rownames(sk.out)<-paste0("D",1:ncol(Phi.est)) 
      }
      rownames(Phi.est)=rownames(Phi.est.unsparse)<-paste0("V",1:p)
      rownames(beta.est)<-colnames(X)
      
      cp.time<-cbind(cp.time,apply(cp.time,1,sum))
      colnames(cp.time)<-c(paste0("D",1:ncol(Phi.est)),"Total")
      
      if(ncol(Phi.est)>1)
      {
        DfD<-diag.level(Y,Phi.est)
        DfD.unsparse<-diag.level(Y,Phi.est.unsparse)
      }else
      {
        DfD<-list(avg.level=1,sub.level=rep(1,n))
        DfD.unsparse<-list(avg.level=1,sub.level=rep(1,n))
      }
    }
    
    if(score.return)
    {
      if(cov.shrinkage)
      {
        colnames(score)=colnames(score.sk)<-paste0("D",1:ncol(Phi.est))
        re<-list(lambda=lambda.out,gamma=Phi.est,beta=beta.est,orthogonality=t(Phi.est)%*%Phi.est,DfD=DfD,shrinkage=sk.out,score=score,score.shrinkage=score.sk,
                 gamma.unsparse=Phi.est.unsparse,DfD.unsparse=DfD.unsparse,orthogonality.unsparse=t(Phi.est.unsparse)%*%Phi.est.unsparse,time=cp.time)
      }else
      {
        colnames(score)<-paste0("D",1:ncol(Phi.est))
        re<-list(lambda=lambda.out,gamma=Phi.est,beta=beta.est,orthogonality=t(Phi.est)%*%Phi.est,DfD=DfD,score=score,
                 gamma.unsparse=Phi.est.unsparse,DfD.unsparse=DfD.unsparse,orthogonality.unsparse=t(Phi.est.unsparse)%*%Phi.est.unsparse,time=cp.time)
      }
    }else
    {
      if(cov.shrinkage)
      {
        re<-list(lambda=lambda.out,gamma=Phi.est,beta=beta.est,orthogonality=t(Phi.est)%*%Phi.est,DfD=DfD,shrinkage=sk.out,
                 gamma.unsparse=Phi.est.unsparse,DfD.unsparse=DfD.unsparse,orthogonality.unsparse=t(Phi.est.unsparse)%*%Phi.est.unsparse,time=cp.time)
      }else
      {
        re<-list(lambda=lambda.out,gamma=Phi.est,beta=beta.est,orthogonality=t(Phi.est)%*%Phi.est,DfD=DfD,
                 gamma.unsparse=Phi.est.unsparse,DfD.unsparse=DfD.unsparse,orthogonality.unsparse=t(Phi.est.unsparse)%*%Phi.est.unsparse,time=cp.time)
      }
    }
    
    return(re)
  }
}
#################################################

###########################################
# Inference of beta
# based on bootstrap
cap_beta_boot<-function(Y,X,gamma=NULL,cov.shrinkage=TRUE,beta=NULL,boot=TRUE,sims=1000,boot.ci.type=c("perc","bca"),conf.level=0.95,verbose=TRUE)
{
  n<-length(Y)
  p<-ncol(Y[[1]])
  Tvec<-rep(NA,n)
  
  q<-ncol(X)
  
  if(is.null(colnames(X)))
  {
    colnames(X)<-c("Intercept",paste0("X",1:(q-1)))
  }
  
  for(i in 1:n)
  {
    Tvec[i]<-nrow(Y[[i]])
  }
  
  if(boot)
  {
    if(is.null(gamma))
    {
      stop("Error! Need gamma value.")
    }else
    {
      beta.boot<-matrix(NA,q,sims)
      
      for(b in 1:sims)
      {
        set.seed(100+b)
        idx.tmp<-sample(1:n,n,replace=TRUE)
        # print(idx.tmp)
        
        Ytmp<-Y[idx.tmp]
        Xtmp<-matrix(X[idx.tmp,],ncol=q)
        
        re.tmp<-NULL
        try(re.tmp<-cap_beta(Ytmp,Xtmp,gamma=gamma,cov.shrinkage=cov.shrinkage,max.itr=max.itr,tol=tol,trace=FALSE,score.return=FALSE))
        if(is.null(re.tmp)==FALSE)
        {
          if(re.tmp$convergence==TRUE)
          {
            beta.boot[,b]<-re.tmp$beta 
          }
        }
        
        if(verbose)
        {
          print(paste0("Bootstrap sample ",b))
        }
      }
      
      for(j in 1:q)
      {
        if(sum(is.na(beta.boot[j,]))>0)
        {
          itmp.nna<-which(is.na(beta.boot[j,])==FALSE)
          dis.cook<-cooks.distance(lm(beta.boot[j,]~1))
          cook.thred<-4/((length(itmp.nna)-2-2))
          itmp<-itmp.nna[which(dis.cook>cook.thred)]
          beta.boot[j,itmp]<-NA
          
          # print(length(itmp))
          
          while(length(itmp)>0)
          {
            itmp.nna<-which(is.na(beta.boot[j,])==FALSE)
            dis.cook<-cooks.distance(lm(beta.boot[j,]~1))
            cook.thred<-4/((length(itmp.nna)-2-2))
            itmp<-itmp.nna[which(dis.cook>cook.thred)]
            beta.boot[j,itmp]<-NA
            
            # print(length(itmp))
          }
        }
      }
      
      beta.est<-apply(beta.boot,1,mean,na.rm=TRUE)
      beta.se<-apply(beta.boot,1,sd,na.rm=TRUE)
      beta.stat<-beta.est/beta.se
      pv<-(1-pnorm(abs(beta.stat)))*2
      
      if(sum(is.na(beta.boot))>0)
      {
        boot.ci.type<-"perc"
      }
      if(boot.ci.type[1]=="bca")
      {
        beta.ci<-t(apply(beta.boot,1,BC.CI,sims=sims,conf.level=conf.level))
      }
      if(boot.ci.type[1]=="perc")
      {
        beta.ci<-t(apply(beta.boot,1,quantile,probs=c((1-conf.level)/2,1-(1-conf.level)/2),na.rm=TRUE))
      }
      
      re<-data.frame(Estiamte=beta.est,SE=beta.se,statistic=beta.stat,pvalue=pv,LB=beta.ci[,1],UB=beta.ci[,2])
      rownames(re)<-colnames(X)
      
      return(list(Inference=re,beta.boot=beta.boot))
    }
  }else
  {
    stop("Error!")
  }
}
###########################################