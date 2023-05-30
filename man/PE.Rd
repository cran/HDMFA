\name{PE}
\alias{PE}
\title{
Projected Estimation for Large-Dimensional Matrix Factor Models 
}
\description{
This function is to fit the matrix factor model via the PE method by projecting the observation matrix onto the row or column factor space.
}
\usage{
PE(X, m1, m2)
}
\arguments{
  \item{X}{
  Input an array with \eqn{T \times p_1 \times p_2}, where \eqn{T} is the sample size, \eqn{p_1} is the the row dimension of each matrix observation and \eqn{p_2} is the the column dimension of each matrix observation.
}
  \item{m1}{
   A positive integer indicating the row factor numbers.
}
  \item{m2}{
  A positive integer indicating the column factor numbers.
}
}
\details{
For the matrix factor models, Yu et al. (2022) propose a projection estimation method to estimate the model parameters. In details, for observations \eqn{\bold{X}_t, t=1,2,\cdots,T}, the data matrix is projected to a lower dimensional space by setting 

\deqn{\bold{Y}_t = \frac{1}{p_2} \bold{X}_t \bold{C}.}

Given \eqn{\bold{Y}_t}, define \deqn{ \bold{M}_1 = \frac{1}{Tp_1} \sum_{t=1}^T \bold{Y}_t \bold{Y}_t^\top, } and then the row factor loading matrix \eqn{\bold{R}} can be estimated by \eqn{\sqrt{p_1}} times the leading \eqn{k_1} eigenvectors of \eqn{\bold{M}_1}. However, the projection matrix \eqn{\bold{C}} is unavailable in practice. A natural solution is to replace it with a consistent initial estimator. The column factor loading matrix \eqn{\bold{C}} can be similarly estimated by projecting \eqn{\bold{X}_t} onto the space of \eqn{\bold{C}} with \eqn{\bold{R}}. See Yu et al. (2022) for the detailed algorithm.
}
\value{
The return value is a list. In this list, it contains the following:
\item{F}{The estimated factor matrix of dimension \eqn{T \times m_1\times m_2}.}
\item{R}{The estimated row loading matrix of dimension \eqn{p_1\times m_1}.}
\item{C}{The estimated column loading matrix of dimension \eqn{p_2\times m_2}.}
}
\references{
Yu, L., He, Y., Kong, X., & Zhang, X. (2022). Projected estimation for large-dimensional matrix factor models. Journal of Econometrics, 229(1), 201-217.
}
\author{
Yong He, Changwei Zhao, Ran Zhao.
}
\examples{
   set.seed(11111)
   T=20;p1=20;p2=20;k1=3;k2=3
   R=matrix(runif(p1*k1,min=-1,max=1),p1,k1)
   C=matrix(runif(p2*k2,min=-1,max=1),p2,k2)
   X=array(0,c(T,p1,p2))
   Y=X;E=Y
   F=array(0,c(T,k1,k2))
   for(t in 1:T){
     F[t,,]=matrix(rnorm(k1*k2),k1,k2)
     E[t,,]=matrix(rnorm(p1*p2),p1,p2)
     Y[t,,]=R\%*\%F[t,,]\%*\%t(C)
   }
   X=Y+E
   
   #Estimate the factor matrices and loadings
   fit=PE(X, k1, k2)
   Rhat=fit$R 
   Chat=fit$C
   Fhat=fit$F
   
   #Estimate the common component
   CC=array(0,c(T,p1,p2))
   for (t in 1:T){
   CC[t,,]=Rhat\%*\%Fhat[t,,]\%*\%t(Chat)
   }
   CC
}

