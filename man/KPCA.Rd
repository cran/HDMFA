\name{KPCA}
\alias{KPCA}
\title{
Estimating the Pair of Factor Numbers via Eigenvalue Ratios Corresponding to \eqn{\alpha}-PCA
}
\description{
The function is to estimate the pair of factor numbers via eigenvalue ratios corresponding to \eqn{\alpha}-PCA. 
}
\usage{
KPCA(X, kmax, alpha = 0)
}
\arguments{
  \item{X}{
   Input an array with \eqn{T \times p_1 \times p_2}, where \eqn{T} is the sample size, \eqn{p_1} is the the row dimension of each matrix observation and \eqn{p_2} is the the column dimension of each matrix observation.
}
  \item{kmax}{
  The user-supplied maximum factor numbers. Here it means the upper bound of the number of row factors and column factors.
}
  \item{alpha}{
   A hyper-parameter balancing the information of the first and second moments    (\eqn{\alpha \geq 1} ). The default is 0.
}
}
\details{
The function \code{\link{KPCA}} uses the eigenvalue-ratio idea to estimate the number of factors. In details, the number of factors \eqn{k_1} is estimated by
\deqn{\hat{k}_1 = \arg \max_{j \leq k_{max}} \frac{\lambda _j (\hat{\bold{M}}_R)}{\lambda _{j+1} (\hat{\bold{M}}_R)},} 
where \eqn{k_{max}} is a given upper bound. \eqn{k_2} is defined similarly with respect to \eqn{\hat{\bold{M}}_C}. See the function \code{\link{alpha_PCA}} for the definition of \eqn{\hat{\bold{M}}_R} and \eqn{\hat{\bold{M}}_C}. For more details, see Chen & Fan (2021).
}
\value{
\item{\eqn{k_1}}{The estimated row factor number.}
\item{\eqn{k_2}}{The estimated column factor number.}
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
   
   KPCA(X, 8, alpha = 0)
}

