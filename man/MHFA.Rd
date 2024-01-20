\name{MHFA}
\alias{MHFA}
\title{
Matrix Huber Factor Analysis
}
\description{
This function is to fit the matrix factor models via the Huber loss. We propose two algorithms to do robust factor analysis. One is based on minimizing the Huber loss of the idiosyncratic error's Frobenius norm, which leads to a weighted iterative projection approach to compute and learn the parameters and thereby named as Robust-Matrix-Factor-Analysis (RMFA). The other one is based on minimizing the element-wise Huber loss, which can be solved by an iterative Huber regression algorithm (IHR).
}
\usage{
MHFA(X, W1=NULL, W2=NULL, m1, m2, method, max_iter = 100, ep = 1e-04)
}
\arguments{
  \item{X}{
  Input an array with \eqn{T \times p_1 \times p_2}, where \eqn{T} is the sample size, \eqn{p_1} is the the row dimension of each matrix observation and \eqn{p_2} is the the column dimension of each matrix observation.
}
  \item{W1}{
  Only if \code{method="E"}, the inital value of row loadings matrix. The default is NULL, which is randomly chosen and all entries from a standard normal distribution.
}
  \item{W2}{
  Only if \code{method="E"}, the inital value of column loadings matrix. The default is NULL, which is randomly chosen and all entries from a standard normal distribution.
}
  \item{m1}{
   A positive integer indicating the row factor numbers.
}
  \item{m2}{
   A positive integer indicating the column factor numbers.
}
  \item{method}{
  Character string, specifying the type of the estimation method to be used.
  \describe{
 \item{\code{"P",}}{indicates minimizing the Huber loss of the idiosyncratic error's Frobenius norm. (RMFA)}
 \item{\code{"E",}}{indicates minimizing the elementwise Huber loss. (IHR)}
 }}
  \item{max_iter}{
  Only if \code{method="E"}, the maximum number of iterations in the iterative Huber regression algorithm. The default is 100.
}
  \item{ep}{
  Only if \code{method="E"}, the stopping critetion parameter in the iterative Huber regression algorithm. The default is \eqn{10^{-4} \times Tp_1 p_2}.
}
}
\details{
For the matrix factor models, He et al. (2021) propose a weighted iterative projection approach to compute and learn the parameters by minimizing the Huber loss function of the idiosyncratic error's Frobenius norm. In details, for observations \eqn{\bold{X}_t, t=1,2,\cdots,T}, define 
\deqn{\bold{M}_c^w = \frac{1}{Tp_2} \sum_{t=1}^T w_t \bold{X}_t \bold{C} \bold{C}^\top \bold{X}_t^\top, \bold{M}_r^w = \frac{1}{Tp_1} \sum_{t=1}^T w_t \bold{X}_t^\top \bold{R} \bold{R}^\top \bold{X}_t.} 
The estimators of loading matrics \eqn{\hat{\bold{R}}} and \eqn{\hat{\bold{C}}} are calculated by \eqn{\sqrt{p_1}} times the leading \eqn{k_1} eigenvectors of \eqn{\bold{M}_c^w} and \eqn{\sqrt{p_2}} times the leading \eqn{k_2} eigenvectors of \eqn{\bold{M}_r^w}.
And \deqn{\hat{\bold{F}}_t=\frac{1}{p_1 p_2}\hat{\bold{R}}^\top \bold{X}_t \hat{\bold{C}}.} For details, see He et al. (2023). 

The other one is based on minimizing the element-wise Huber loss. Define 
\deqn{M_{i,Tp_2}(\bold{r}, \bold{F}_t, \bold{C})=\frac{1}{Tp_2} \sum_{t=1}^{T} \sum_{j=1}^{p_2} H_\tau \left(x_{t,ij}-\bold{r}_i^\top\bold{F}_t\bold{c}_j \right),}
\deqn{M_{i,Tp_1}(\bold{R}, \bold{F}_t, \bold{c})=\frac{1}{Tp_1}\sum_{t=1}^T\sum_{i=1}^{p_1} H_\tau \left(x_{t,ij}-\bold{r}_i^\top\bold{F}_t\bold{c}_j\right),}
\deqn{M_{t,p_1 p_2}(\bold{R}, \mathrm{vec}(\bold{F}), \bold{C})=\frac{1}{p_1 p_2} \sum_{i=1}^{p_1}\sum_{j=1}^{p_2} H_\tau \left(x_{t,ij}-(\bold{c}_j \otimes \bold{r}_i)^\top \mathrm{vec}(\bold{F})\right).} This can be seen as Huber regression as each time optimizing one argument while keeping the other two fixed.
}
\value{
The return value is a list. In this list, it contains the following:
\item{F}{The estimated factor matrix of dimension \eqn{T \times m_1\times m_2}.}
\item{R}{The estimated row loading matrix of dimension \eqn{p_1\times m_1}, satisfying \eqn{\bold{R}^\top\bold{R}=p_1\bold{I}_{m_1}}.}
\item{C}{The estimated column loading matrix of dimension \eqn{p_2\times m_2}, satisfying \eqn{\bold{C}^\top\bold{C}=p_2\bold{I}_{m_2}}.}
}
\references{
He, Y., Kong, X., Yu, L., Zhang, X., & Zhao, C. (2023). Matrix factor analysis: From least squares to iterative projection. Journal of Business & Economic Statistics, 1-26.

He, Y., Kong, X. B., Liu, D., & Zhao, R. (2023). Robust Statistical Inference for Large-dimensional Matrix-valued Time Series via Iterative Huber Regression. <arXiv:2306.03317>.
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
   
   #Estimate the factor matrices and loadings by RMFA
   fit1=MHFA(X, m1=3, m2=3, method="P")
   Rhat1=fit1$R 
   Chat1=fit1$C
   Fhat1=fit1$F
   
   #Estimate the factor matrices and loadings by IHR
   fit2=MHFA(X, W1=NULL, W2=NULL, 3, 3, "E")
   Rhat2=fit2$R 
   Chat2=fit2$C
   Fhat2=fit2$F
   
   #Estimate the common component by RMFA
   CC1=array(0,c(T,p1,p2))
   for (t in 1:T){
   CC1[t,,]=Rhat1\%*\%Fhat1[t,,]\%*\%t(Chat1)
   }
   CC1
   
   #Estimate the common component by IHR
   CC2=array(0,c(T,p1,p2))
   for (t in 1:T){
   CC2[t,,]=Rhat2\%*\%Fhat2[t,,]\%*\%t(Chat2)
   }
   CC2
}
