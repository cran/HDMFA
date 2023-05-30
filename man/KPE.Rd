\name{KPE}
\alias{KPE}
\title{
Estimating the Pair of Factor Numbers via Eigenvalue Ratios Corresponding to PE
}
\description{
The function is to estimate the pair of factor numbers via eigenvalue ratios corresponding to PE method. 
}
\usage{
KPE(X, kmax, c = 0)
}
\arguments{
  \item{X}{
  Input an array with \eqn{T \times p_1 \times p_2}, where \eqn{T} is the sample size, \eqn{p_1} is the the row dimension of each matrix observation and \eqn{p_2} is the the column dimension of each matrix observation.
}
  \item{kmax}{
  The user-supplied maximum factor numbers. Here it means the upper bound of the number of row factors and column factors.
}
  \item{c}{
  A constant to avoid vanishing denominators. The default is 0.
}
}
\details{
The function \code{\link{KPE}} uses the eigenvalue-ratio idea to estimate the number of factors. 
First, obtain the initial estimators \eqn{\hat{\bold{R}}} and \eqn{\hat{\bold{C}}}. Second, define 
\deqn{\hat{\bold{Y}}_t=\frac{1}{p_2}\bold{X}_t\hat{\bold{C}}, \hat{\bold{Z}}_t=\frac{1}{p_1}\bold{X}_t^\top\hat{\bold{R}},} and 
\deqn{\tilde{\bold{M}}_1=\frac{1}{Tp_1}\hat{\bold{Y}}_t\hat{\bold{Y}}_t^\top, \tilde{\bold{M}}_2=\frac{1}{Tp_2}\sum_{t=1}^T\hat{\bold{Z}}_t\hat{\bold{Z}}_t^\top,} the number of factors \eqn{k_1} is estimated by
\deqn{\hat{k}_1 = \arg \max_{j \leq k_{max}} \frac{\lambda_j (\tilde{\bold{M}}_1)}{\lambda _{j+1} (\tilde{\bold{M}}_1)},}
where \eqn{k_{max}} is a predetermined upper bound for \eqn{k_1}. The estimation of \eqn{k_2} is defined similarly with respect to \eqn{\tilde{\bold{M}}_2}. 
For details, see Yu et al. (2022).
}
\value{
\item{\eqn{k_1}}{The estimated row factor number.}
\item{\eqn{k_2}}{The estimated column factor number.}
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
   
   KPE(X, 8, c = 0)
}