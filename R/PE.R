PE <- function(X,m1,m2){
  T=dim(X)[1]
  p1=dim(X)[2]
  p2=dim(X)[3]
  
  C0 <- alpha_PCA(X,m1,m2,alpha=0)$C
  R0 <- alpha_PCA(X,m1,m2,alpha=0)$R
  Rhat=SCR2(X,m1,C0)*sqrt(p1)
  Chat=SCC2(X,m2,R0)*sqrt(p2)
  Fhat=array(0,c(T,m1,m2))
  for(t in 1:T){
    Fhat[t,,]=t(Rhat)%*%X[t,,]%*%Chat/(p1*p2)
  }
  return(list(R=Rhat,C=Chat,F=Fhat))
}

KPE <- function(X,kmax,c=0){
  T=dim(X)[1]
  p1=dim(X)[2]
  p2=dim(X)[3]
  d1=(1/sqrt(T*p1)+1/sqrt(T*p2)+1/p2)*c
  d2=(1/sqrt(T*p1)+1/sqrt(T*p2)+1/p1)*c
  
  M1=matrix(0,p1,p1)
  for(t in 1:T){
    M1=M1+X[t,,]%*%t(X[t,,])
  }
  R1=svds(M1,kmax,kmax,kmax)$u
  M2=matrix(0,p2,p2)
  for(t in 1:T){
    M2=M2+t(X[t,,])%*%X[t,,]
  }
  C1=svds(M2,kmax,kmax,kmax)$u
  
  k2=0
  k1=0
  iter=0
  k1_new=kmax
  k2_new=kmax
  while ((k1_new != k1 |k2_new!=k2) & iter<10) {
    k1=k1_new
    k2=k2_new
    M2=matrix(0,p2,p2)
    for(t in 1:T){
      M2=M2+t(X[t,,])%*%R1[,1:k1]%*%t(R1[,1:k1])%*%X[t,,]
    }
    eigval=eigen(M2)$values[1:(1+kmax)]
    k2_new=which.max(eigval[1:kmax]/(eigval[2:(1+kmax)]+d2))
    
    M1=matrix(0,p1,p1)
    for(t in 1:T){
      M1=M1+X[t,,]%*%C1[,1:k2_new]%*%t(C1[,1:k2_new])%*%t(X[t,,])
    }
    eigval=eigen(M1)$values[1:(1+kmax)]
    k1_new=which.max(eigval[1:kmax]/(eigval[2:(1+kmax)]+d1))
    iter=iter+1
  }
  return(c(k1,k2))
}
