\name{alpha_PCA}
\alias{alpha_PCA}
\title{
Statistical Inference for High-Dimensional Matrix-Variate Factor Model
}
\description{
This function is to fit the matrix factor model via the \eqn{\alpha}-PCA method by conducting eigen-analysis of a weighted average of the sample mean and the column (row) sample covariance matrix through a hyper-parameter \eqn{\alpha}.
}
\usage{
alpha_PCA(X, m1, m2, alpha = 0)
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
  \item{alpha}{
  A hyper-parameter balancing the information of the first and second moments    (\eqn{\alpha \geq -1} ). The default is 0.
}
}
\details{
For the matrix factor models, Chen & Fan (2021) propose an estimation procedure, i.e. \eqn{\alpha}-PCA. The method aggregates the information in both first and second moments and extract it via a spectral method. In detail, for observations \eqn{\bold{X}_t, t=1,2,\cdots,T}, define 
\deqn{\hat{\bold{M}}_R = \frac{1}{p_1 p_2} \left( (1+\alpha) \bar{\bold{X}} \bar{\bold{X}}^\top + \frac{1}{T} \sum_{t=1}^T (\bold{X}_t - \bar{\bold{X}}) (\bold{X}_t - \bar{\bold{X}})^\top \right),}
\deqn{\hat{\bold{M}}_C = \frac{1}{p_1 p_2} \left( (1+\alpha) \bar{\bold{X}}^\top \bar{\bold{X}} + \frac{1}{T} \sum_{t=1}^T (\bold{X}_t - \bar{\bold{X}})^\top (\bold{X}_t - \bar{\bold{X}}) \right),} 
where \eqn{\alpha \in} [-1,\eqn{+\infty}], \eqn{\bar{\bold{X}} = \frac{1}{T} \sum_{t=1}^T \bold{X}_t}, \eqn{\frac{1}{T} \sum_{t=1}^T (\bold{X}_t - \bar{\bold{X}}) (\bold{X}_t - \bar{\bold{X}})^\top} and \eqn{\frac{1}{T} \sum_{t=1}^T (\bold{X}_t - \bar{\bold{X}})^\top (\bold{X}_t - \bar{\bold{X}})} are the sample row and column covariance matrix, respectively. The loading matrices \eqn{\bold{R}} and \eqn{\bold{C}} are estimated as \eqn{\sqrt{p_1}} times the top \eqn{k_1} eigenvectors of \eqn{\hat{\bold{M}}_R} and \eqn{\sqrt{p_2}} times the top \eqn{k_2} eigenvectors of \eqn{\hat{\bold{M}}_C}, respectively. For details, see Chen & Fan (2021). 
}
\value{
The return value is a list. In this list, it contains the following:
\item{F}{The estimated factor matrix of dimension \eqn{T \times m_1\times m_2}.}
\item{R}{The estimated row loading matrix of dimension \eqn{p_1\times m_1}, satisfying \eqn{\bold{R}^\top\bold{R}=p_1\bold{I}_{m_1}}.}
\item{C}{The estimated column loading matrix of dimension \eqn{p_2\times m_2}, satisfying \eqn{\bold{C}^\top\bold{C}=p_2\bold{I}_{m_2}}.}
}
\references{
Chen, E. Y., & Fan, J. (2021). Statistical inference for high-dimensional matrix-variate factor models. Journal of the American Statistical Association, 1-18.
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
   fit=alpha_PCA(X, k1, k2, alpha = 0)
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
