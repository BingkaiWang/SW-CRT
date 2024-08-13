########################################
# continuous outcome GEE
# Bias-corrected Variance
########################################

########################################
# Input
# y: N by 1 vector of outcomes
# X: N by p design matrix (including intercept)
# beta: p by 1 mean model parameter estimates
# tau: 1 by 1 within-period correlation parameter estimate
# rho: 1 by 1 decay correlation parameter estimate
# phi: 1 by 1 dispersion parameter estimate
# id: N by 1 cluster indicator
# t: period
# mv: vector of cluster-period sizes
########################################

contGEEBCV=function(y,X,beta,tau,rho,phi,id,t,mv){
  
  # Creates two vectors that have the start and end points for each cluster
  BEGINEND=function(n){
    last=cumsum(n)
    first=last-n+1
    return(cbind(first,last))
  }
  
  # Calculate the inverse for symmetric and positive finite matrix
  FINDINV=function(A){
    require(MASS)
    AHALF=chol(A)
    GINV=ginv(AHALF)
    AINV=tcrossprod(GINV)
    return(AINV)
  }
  
  # Generate exponential decay correlation
  ed <- function(ni, t, tau, rho){
    Ni <- sum(ni)
    R<-matrix(0,Ni,Ni)
    id<-BEGINEND(ni)
    for(j in 1:(t-1)){
      row=col=(id[j,1]:id[j,2])
      R[row,col]=(1-tau)*diag(ni[j])+tau*matrix(1,ni[j],ni[j])
      for(k in (j+1):t){
        col=(id[k,1]:id[k,2])
        R[row,col]=tau * rho^(abs(k-j))
      }
    }
    row=col=(id[t,1]:id[t,2])
    R[row,col]=(1-tau)*diag(ni[t])+tau*matrix(1,ni[t],ni[t])
    R[lower.tri(R)] = t(R)[lower.tri(R)]
    
    return(R)
  }
  
  # Score function
  SCORE=function(beta,tau,rho,phi,y,X,n,p,t,mv){
    U=rep(0,p)
    UUtran=Ustar=matrix(0,p,p)
    locx=BEGINEND(n)
    
    for(i in 1:length(n)){
      X_c=X[locx[i,1]:locx[i,2],,drop=FALSE]
      y_c=y[locx[i,1]:locx[i,2]]
      
      U_c=rep(0,p)
      Ustar_c=matrix(0,p,p)
      mu_c=X_c%*%beta
      
      C=X_c
      A=y_c-mu_c
      
      mviv <- mv[(t*i-(t-1)):(t*i)]
      R <- ed(mviv, t, tau, rho)
      INVR=solve(R)
      INVB=INVR/phi
      
      U_c=t(C)%*%INVB%*%A
      UUtran_c=tcrossprod(U_c)
      Ustar_c=t(C)%*%INVB%*%C
      U=U+U_c
      UUtran=UUtran+UUtran_c
      Ustar=Ustar+Ustar_c
    }
    return(list(U=U,UUtran=UUtran,Ustar=Ustar))
  }
  
  # Compute (A - mm`)^{-1}c without performing the inverse directly
  INVBIG=function(ainvc,ainvm,m,c,start,end){
    for(i in start:end){
      b=ainvm[,i]
      bt=t(b)
      btm=bt%*%m
      btmi=btm[,i]
      gam=1-btmi
      bg=b/gam
      ainvc=ainvc+bg%*%(bt%*%c)
      if(i<end){
        ainvm=ainvm+bg%*%btm
      }
    }
    return(ainvc)
  }
  
  # Creates bias-corrected covariance matrix of beta
  p=ncol(X)
  n=as.numeric(table(id))
  SCORE_RES=SCORE(beta,tau,rho,phi,y,X,n,p,t,mv)
  U=SCORE_RES$U
  UUtran=SCORE_RES$UUtran
  Ustar=SCORE_RES$Ustar
  
  # Naive or Model-based estimator
  naive=FINDINV(Ustar)
  
  # BC0 or usual Sandwich estimator     
  robust=naive%*%UUtran%*%t(naive)
  
  # new commands to compute INV(I - H1)
  eigenRES1=eigen(naive)
  evals1=eigenRES1$values
  evecs1=eigenRES1$vectors
  sqrevals1=sqrt(evals1)
  sqe1=evecs1%*%diag(sqrevals1)
  
  # Bias-corrected variance
  Ustar_c_array=UUtran_c_array=array(0,c(p,p,length(n)))
  UUtran=UUbc=UUbc2=UUbc3=Ustar=matrix(0,p,p)
  
  locx=BEGINEND(n)
  
  for(i in 1:length(n)){
    X_c=X[locx[i,1]:locx[i,2],,drop=FALSE]
    y_c=y[locx[i,1]:locx[i,2]]
    mu_c=X_c%*%beta
    
    U_i=U_c=rep(0,p)
    Ustar_c=matrix(0,p,p)
    
    # commands for beta
    C=X_c
    A=y_c-mu_c
    
    mviv <- mv[(t*i-(t-1)):(t*i)]
    R <- ed(mviv, t, tau, rho)
    INVR=solve(R)
    INVB=INVR/phi
    U_i=t(C)%*%INVB%*%A
    
    # commands for generalized inverse - beta
    ai1=INVB
    mm1=C%*%sqe1
    ai1A=ai1%*%A
    ai1m1=ai1%*%mm1
    ai1A=INVBIG(ai1A,ai1m1,mm1,A,1,p)
    U_c=t(C)%*%ai1A
    
    Ustar_c=t(C)%*%INVB%*%C
    Ustar=Ustar+Ustar_c
    UUtran_c=tcrossprod(U_i)
    UUtran=UUtran+UUtran_c
    UUbc_c=tcrossprod(U_c)
    UUbc=UUbc+UUbc_c
    UUbc_ic=tcrossprod(U_c,U_i)
    UUbc2=UUbc2+UUbc_ic
    
    Ustar_c_array[,,i]=Ustar_c
    UUtran_c_array[,,i]=UUtran_c
  }
  
  # calculating adjustment factor for BC3
  # for(i in 1:length(n)){      
  #   Hi=diag(1/sqrt(1-pmin(0.75,c(diag(Ustar_c_array[,,i]%*%naive)))))
  #   UUbc3=UUbc3+Hi%*%UUtran_c_array[,,i]%*%Hi
  # }
  
  # BC1 or Variance estimator due to Kauermann and Carroll (2001);
  varKC=naive%*%(UUbc2+t(UUbc2))%*%t(naive)/2
  
  # BC2 or Variance estimator due to Mancl and DeRouen (2001);
  varMD=naive%*%UUbc%*%t(naive)
  
  # BC3 or Variance estimator due to Fay and Graubard (2001);
  # varFG=naive%*%UUbc3%*%t(naive)
  
  ########################################
  # Output
  # naive: naive or model-based var
  # robust: robust sandwich var
  # varMD: bias-corrected sandwich var due to Mancl and DeRouen (2001)
  # varKC: bias-corrected sandwich var due to Kauermann and Carroll (2001)
  # varFG: bias-corrected sandwich var due to Fay and Graubard (2001)
  ########################################
  return(list(naive=naive,robust=robust,varMD=varMD,varKC=varKC))
}