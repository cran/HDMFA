RMFA <- function(X,k1,k2){
  T=dim(X)[1]
  p1=dim(X)[2]
  p2=dim(X)[3]
  
  R0=alpha_PCA(X,k1,k2,alpha=0)$R
  C0=alpha_PCA(X,k1,k2,alpha=0)$C
  
  tau=c()
  for (i in 1:T) {
    tau[i]=sqrt(abs(sum(diag(t(X[i,,])%*%X[i,,]-
                               t(X[i,,])%*%R0%*%t(R0)%*%X[i,,]%*%C0%*%t(C0)/(p1*p2)))))
  }
  tau=median(tau)
  
  w=c();r=c()
  for (i in 1:T) {
    r[i]=sqrt(abs(sum(diag(t(X[i,,])%*%X[i,,]-
                             t(X[i,,])%*%R0%*%t(R0)%*%X[i,,]%*%C0%*%t(C0)/(p1*p2)))))
    if(r[i]<=tau){
      w[i]=1/2
    }else{
      w[i]=tau/2/r[i]
    }
  }
  R=R0;C=C0
  
  L_init=Huberloss(X,R,C,tau)
  L_old=L_init+1
  L_new=L_init
  iter=0
  
  while (L_old>L_new) {
    Fr=rowcov(X,C,w)*T*p2
    Fr_eigen=eigen(Fr)
    R_new=Fr_eigen$vectors[,1:k1]
    Fc=colcov(X,R,w)*T*p1
    Fc_eigen=eigen(Fc)
    C_new=Fc_eigen$vectors[,1:k2]
    tau=c()
    for (i in 1:T) {
      tau[i]=sqrt(abs(sum(diag(t(X[i,,])%*%X[i,,]-
                                 t(X[i,,])%*%R_new%*%t(R_new)%*%X[i,,]%*%C_new%*%t(C_new)/(p1*p2)))))
    }
    tau=median(tau)
    L_old=L_new
    L_new=Huberloss(X,R_new,C_new,tau)
    #weights
    w=c();r=c()
    for (i in 1:T) {
      r[i]=sqrt(abs(sum(diag(t(X[i,,])%*%X[i,,]-
                               t(X[i,,])%*%R_new%*%t(R_new)%*%X[i,,]%*%C_new%*%t(C_new)/(p1*p2)))))
      if(r[i]<=tau){w[i]=1/2}else{w[i]=tau/2/r[i]}
    }
  }
  Fhat=array(0,c(T,k1,k2))
  for (t in 1:T) {
    Fhat[t,,]=t(R_new)%*%X[t,,]%*%C_new/(p1*p2)
  }
  
  return(list(R=R_new,C=C_new,F=Fhat))
}

KRMFM <- function(X,kmax,c=0){
  T=dim(X)[1]
  p1=dim(X)[2]
  p2=dim(X)[3]
  d1=(1/sqrt(T*p1)+1/sqrt(T*p2)+1/p2)*c
  d2=(1/sqrt(T*p1)+1/sqrt(T*p2)+1/p1)*c
  
  A=matrix(0,p1,p1)
  for(t in 1:T){
    A=A+X[t,,]%*%t(X[t,,])
  }
  R1=svds(A,kmax,kmax,kmax)$u
  B=matrix(0,p2,p2)
  for(t in 1:T){
    B=B+t(X[t,,])%*%X[t,,]
  }
  C1=svds(B,kmax,kmax,kmax)$u
  
  tau=c()
  for (i in 1:T) {
    tau[i]=sqrt(abs(sum(diag(t(X[i,,])%*%X[i,,]-
                               t(X[i,,])%*%R1%*%t(R1)%*%X[i,,]%*%C1%*%t(C1)/(p1*p2)))))
  }
  tau=median(tau)
  #######weights
  w=c();r=c()
  for (i in 1:T) {
    r[i]=sqrt(abs(sum(diag(t(X[i,,])%*%X[i,,]-
                             t(X[i,,])%*%R1%*%t(R1)%*%X[i,,]%*%C1%*%t(C1)/(p1*p2)))))
    if(r[i]<=tau){
      w[i]=1/2
    }else{
      w[i]=tau/2/r[i]
    }
  }
  
  k2=0;k1=0
  iter=0
  k1_new=kmax;k2_new=kmax
  while ((k1_new != k1 |k2_new!=k2) & iter<10) {
    k1=k1_new
    k2=k2_new
    B=matrix(0,p2,p2)
    for(t in 1:T){
      B=B+w[t]*(t(X[t,,])%*%R1[,1:k1]%*%t(R1[,1:k1])%*%X[t,,])
    }
    eigval=eigen(B)$values[1:(1+kmax)]
    k2_new=which.max(eigval[1:kmax]/(eigval[2:(1+kmax)]+d2))
    
    A=matrix(0,p1,p1)
    for(t in 1:T){
      A=A+w[i]*(X[t,,]%*%C1[,1:k2_new]%*%t(C1[,1:k2_new])%*%t(X[t,,]))
    }
    eigval=eigen(A)$values[1:(1+kmax)]
    k1_new=which.max(eigval[1:kmax]/(eigval[2:(1+kmax)]+d1))
    iter=iter+1
  }
  return(c(k1,k2))
}

IHR <- function(X,W1,W2,m1,m2,max_iter=100,ep=1e-4){
  T=dim(X)[1]
  p1=dim(X)[2]
  p2=dim(X)[3]
  ep=ep*T*p1*p2
  
  if(class(W1)[1]=="NULL"&class(W2)[1]=="NULL"){
    W1=matrix(rnorm(p1*m1),p1,m1); W1=sqrt(p1)*svd(W1)$u
    W2=matrix(rnorm(p2*m2),p2,m2); W2=sqrt(p2)*svd(W2)$u
  }else{
    W1=W1;W2=W2
  }
  #The initial estimate of factor matrix
  Fhat=array(0,c(T,m1,m2)) 
  for(t in 1:T){
    Fhat[t,,]=t(W1)%*%X[t,,]%*%W2/(p1*p2)
  }
  C_old=W2
  R_old=W1
  Fhat_old=Fhat
  
  for (m in 1:max_iter) {
    #solve the M_{i,Tp2} to get Rhat
    A_r=c()
    for (j in 1:p2) {
      for (t in 1:T) {
        A_r=rbind(A_r,t(C_old[j,])%*%t(matrix(Fhat_old[t,,],m1,m2)))
      }
    }
    Z_r=matrix(0,T*p2,p1)
    R_tilde=matrix(0,p1,m1)
    for (i in 1:p1) {
      Z_r[,i]=c(X[,i,])
      R_tilde[i,]=rlm(y=Z_r[,i],x=A_r,method="M",scale.est="Huber")$coef
    }
    
    #Solve the M_{i,Tp1} to get Chat
    A_c=c()
    for (i in 1:p1) {
      for (t in 1:T) {
        A_c=rbind(A_c,R_tilde[i,]%*%Fhat_old[t,,])
      }
    }
    Z_c=matrix(0,T*p1,p2)
    C_tilde=matrix(0,p2,m2)
    for (j in 1:p2) {
      Z_c[,j]=c(X[,,j])
      C_tilde[j,]=rlm(y=Z_c[,j],x=A_c,method="M",scale.est="Huber")$coef
    }
    
    #solve M_{t,p1p2} to get Fhat
    A_f=c()
    for (j in 1:p2) {
      for (i in 1:p1) {
        A_f=rbind(A_f,kronecker(C_tilde[j,],R_tilde[i,]))
      }
    }
    Fhat_new1=array(0,c(T,m1,m2))
    Z_f=matrix(0,p1*p2,T)
    F_tilde=matrix(0,T,m1*m2)
    for (t in 1:T) {
      Z_f[,t]=c(X[t,,])
      F_tilde[t,]=rlm(y=Z_f[,t],x=A_f,method="M",scale.est="Huber")$coef
      Fhat_new1[t,,]=matrix(F_tilde[t,],m1,m2,byrow=FALSE)
    }
    
    R_new=Normalization(X,R_tilde,C_tilde,Fhat_new1,m1,m2)$R
    C_new=Normalization(X,R_tilde,C_tilde,Fhat_new1,m1,m2)$C
    Fhat_new=Normalization(X,R_tilde,C_tilde,Fhat_new1,m1,m2)$F
    
    CC_old=array(0,c(T,p1,p2));CC_new=array(0,c(T,p1,p2))
    CCdiff=c()
    for(t in 1:T){
      CC_new[t,,]=R_new%*%Fhat_new[t,,]%*%t(C_new)
      CC_old[t,,]=R_old%*%Fhat_old[t,,]%*%t(C_old)
      CCdiff[t]=norm(CC_new[t,,]-CC_old[t,,],"F")
    }
    
    if(sum(CCdiff) <= ep){
      return(list(R=R_new,C=C_new,F=Fhat_new,iter=m))
    }else{
      C_old=C_new
      R_old=R_new
      Fhat_old=Fhat_new
    }
  }
  return(list(R=R_new,C=C_new,F=Fhat_new,iter=max_iter))
}

KIHR <- function(X,W1,W2,kmax,method,max_iter=100,c=1e-4,ep=1e-4){
  T=dim(X)[1];p1=dim(X)[2];p2=dim(X)[3]
  
  D=min(sqrt(T*p1),sqrt(T*p2),sqrt(p1*p2))
  Fhat=IHR(X,W1,W2,kmax,kmax,max_iter,ep=1e-4)$F
  sigma1=matrix(0,kmax,kmax);sigma2=matrix(0,kmax,kmax)
  for (t in 1:T) {
    sigma1=sigma1+Fhat[t,,]%*%t(Fhat[t,,])/T
    sigma2=sigma2+t(Fhat[t,,])%*%Fhat[t,,]/T
  }
  if (method=="E_RM"){
    k1=sum(diag(sigma1)>sigma1[1,1]*D^(-2/3))
    k2=sum(diag(sigma2)>sigma2[1,1]*D^(-2/3))
  }
  if (method=="E_ER"){
    eigval1=diag(sigma1)
    k1=which.max(eigval1[1:(kmax-1)]/(eigval1[2:kmax]+c*D^(-2)))
    eigval2=diag(sigma2)
    k2=which.max(eigval2[1:(kmax-1)]/(eigval2[2:kmax]+c*D^(-2)))
  }
  return(c(k1,k2))
}

MHFA <- function(X,W1=NULL,W2=NULL,m1,m2,method,max_iter=100,ep=1e-4){
  if (method=="P"){
    fit=RMFA(X,m1,m2)
    R=fit$R
    C=fit$C
    F=fit$F
    return(list(R=R,C=C,F=F))
  }
  if (method=="E"){
    fit=IHR(X,W1,W2,m1,m2,max_iter=100,ep=1e-4)
    R=fit$R
    C=fit$C
    F=fit$F
    return(list(R=R,C=C,F=F))
  }
}

KMHFA <- function(X,W1=NULL,W2=NULL,kmax,method,max_iter=100,c=1e-4,ep=1e-4){
  if (method=="P"){
    K=KRMFM(X,kmax,c)
    k1=K[1];k2=K[2]
    return(list(k1=k1,k2=k2))
  }
  if (method=="E_RM"){
    K=KIHR(X,W1,W2,kmax,"E_RM",max_iter,c,ep)
    k1=K[1];k2=K[2]
    return(list(k1=k1,k2=k2))
  }
  if (method=="E_ER"){
    K=KIHR(X,W1,W2,kmax,"E_ER",max_iter,c,ep)
    k1=K[1];k2=K[2]
    return(list(k1=k1,k2=k2))
  }
}
