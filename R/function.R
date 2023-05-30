#Calculate the weighted row-row and col-col covariance matrices
rowcov <- function(X, C, Weights) {
  T=dim(X)[1]
  p1=dim(X)[2]
  Fr=matrix(0,p1,p1)
  C2=C%*%t(C)
  for (t in 1:T) {
    Fr=Fr+Weights[t]*(X[t,,]%*%C2%*%t(X[t,,]))
  }
  return(Fr)
}

colcov <- function(X, R, Weights){
  T=dim(X)[1]
  p2=dim(X)[3]
  Fc=matrix(0,p2,p2)
  R2=R%*%t(R)
  for (i in 1:T) {
    Fc=Fc+Weights[i]*(t(X[i,,])%*%R2%*%X[i,,])
  }
  return(Fc)
}

#Calculate Huber loss
Huberloss <- function(X,R,C,tau){
  T=dim(X)[1];p1=dim(X)[2];p2=dim(X)[3]
  r=c();loss=c()
  T=dim(X)[1]
  R2=R%*%t(R)
  C2=C%*%t(C)
  for (i in 1:T) {
    r[i]=sqrt(abs(sum(diag(t(X[i,,])%*%X[i,,]-
                             t(X[i,,])%*%R2%*%X[i,,]%*%C2/(p1*p2)))))
    if(r[i]<=tau){
      loss[i]=(r[i]^2)/2
    }else{
      loss[i]=tau*abs(r[i])-tau^2/2
    }
  }
  result=sum(loss)
  return(result)
}

#alpha PCA，W1 W2 are the first and second moment terms in their method
aPCAR<-function(W1,W2,k1,alpha){
  M=alpha*W1+W2
  return(svds(M,k1,k1,k1)$u)
}

aPCAC<-function(S1,S2,k2,alpha){
  M=alpha*S1+S2
  return(svds(M,k2,k2,k2)$u)
}

#Sample covariance with 2 step
SCR2<-function(Y,k1,C0){
  T=dim(Y)[1]
  p1=dim(Y)[2]
  p2=dim(Y)[3]
  M=matrix(0,p1,p1)
  for(t in 1:T){
    M=M+Y[t,,]%*%C0%*%t(C0)%*%t(Y[t,,])
  }
  return(svds(M,k1,k1,k1)$u)
}

SCC2<-function(Y,k2,R0){
  T=dim(Y)[1]
  p1=dim(Y)[2]
  p2=dim(Y)[3]
  M=matrix(0,p2,p2)
  for(t in 1:T){
    M=M+t(Y[t,,])%*%R0%*%t(R0)%*%Y[t,,]
  }
  return(svds(M,k2,k2,k2)$u)
}

Normalization <- function(X,R0,C0,F0,m1,m2){
  T=dim(X)[[1]]
  p1=dim(X)[[2]]
  p2=dim(X)[[3]]
  
  U_r=svd(R0)$u;svd_R0=svd(R0)$d
  Q_r=diag(svd_R0,length(svd_R0))%*%t(svd(R0)$v)
  U_c=svd(C0)$u;svd_C0=svd(C0)$d
  Q_c=diag(svd_C0,length(svd_C0))%*%t(svd(C0)$v)
  Sigma1=matrix(0,m1,m1);Sigma2=matrix(0,m2,m2)
  for (t in 1:T) {
    Sigma1=Sigma1+Q_r%*%F0[t,,]%*%t(C0)%*%C0%*%t(matrix(F0[t,,],m1,m2))%*%t(Q_r)/(T*p1*p2)
    Sigma2=Sigma2+Q_c%*%t(matrix(F0[t,,],m1,m2))%*%t(R0)%*%R0%*%F0[t,,]%*%t(Q_c)/(T*p1*p2)
  }
  Gamma1=eigen(Sigma1)$vectors
  Gamma2=eigen(Sigma2)$vectors
  
  Rnew=sqrt(p1)*U_r%*%Gamma1
  Cnew=sqrt(p2)*U_c%*%Gamma2
  Fnew=array(0,c(T,m1,m2))
  for (t in 1:T) {
    Fnew[t,,]=t(Gamma1)%*%Q_r%*%F0[t,,]%*%t(Q_c)%*%Gamma2/(sqrt(p1*p2))
  }
  return(list(R=Rnew,C=Cnew,F=Fnew))
}
