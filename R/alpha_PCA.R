alpha_PCA <- function(X,m1,m2,alpha=0){
  T=dim(X)[1]
  p1=dim(X)[2]
  p2=dim(X)[3]
  
  barX=apply(X, c(2,3), mean)
  MR1=barX%*%t(barX)
  MR2=matrix(0,p1,p1)
  for(t in 1:T){
    MR2=MR2+X[t,,]%*%t(X[t,,])
  }
  MR2=MR2/T
  
  MC1=t(barX)%*%barX
  MC2=matrix(0,p2,p2)
  for(t in 1:T){
    MC2=MC2+t(X[t,,])%*%X[t,,]
  }
  MC2=MC2/T
  
  Rhat=aPCAR(MR1,MR2,m1,alpha)
  Chat=aPCAC(MC1,MC2,m2,alpha)
  Fhat=array(0,c(T,m1,m2))
  for(t in 1:T){
    Fhat[t,,]=t(Rhat)%*%X[t,,]%*%Chat/(p1*p2)
  }
  return(list(R=Rhat,C=Chat,F=Fhat))
}

KPCA <- function(X,kmax,alpha=0){
  T=dim(X)[1]
  p1=dim(X)[2]
  p2=dim(X)[3]
  
  barX=apply(X,c(2,3),mean)
  W1=barX%*%t(barX)
  W2=matrix(0,p1,p1)
  for(t in 1:T){
    W2=W2+X[t,,]%*%t(X[t,,])
  }
  W2=W2/T
  
  S1=t(barX)%*%barX
  S2=matrix(0,p2,p2)
  for(t in 1:T){
    S2=S2+t(X[t,,])%*%X[t,,]
  }
  S2=S2/T
  M1=alpha*W1+W2
  v1=svds(M1,kmax+1,kmax+1,kmax+1)$d
  k1=which.max(v1[1:kmax]/v1[2:(kmax+1)])
  
  M2=alpha*S1+S2
  v2=svds(M2,kmax+1,kmax+1,kmax+1)$d
  k2=which.max(v2[1:kmax]/v2[2:(kmax+1)])
  return(c(k1,k2))
}