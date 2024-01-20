\name{KMHFA}
\alias{KMHFA}
\title{
Estimating the Pair of Factor Numbers via Eigenvalue Ratios or Rank Minimization.
}
\description{
The function is to estimate the pair of factor numbers via eigenvalue-ratio corresponding to RMFA method or rank minimization and eigenvalue-ratio corresponding to Iterative Huber Regression (IHR). 
}
\usage{
KMHFA(X, W1 = NULL, W2 = NULL, kmax, method, max_iter = 100, c = 1e-04, ep = 1e-04)
}
\arguments{
  \item{X}{
  Input an array with \eqn{T \times p_1 \times p_2}, where \eqn{T} is the sample size, \eqn{p_1} is the the row dimension of each matrix observation and \eqn{p_2} is the the column dimension of each matrix observation.
}
  \item{W1}{
  Only if \code{method="E_RM"} or \code{method="E_ER"}, the inital value of row loadings matrix. The default is NULL, which is randomly chosen and all entries from a standard normal distribution.
}
  \item{W2}{
  Only if \code{method="E_RM"} or \code{method="E_ER"}, the inital value of column loadings matrix. The default is NULL, which is randomly chosen and all entries from a standard normal distribution.
}
  \item{kmax}{
  The user-supplied maximum factor numbers. Here it means the upper bound of the number of row factors and column factors.
}
  \item{method}{
   Character string, specifying the type of the estimation method to be used.
  \describe{
  \item{\code{"P",}}{the robust iterative eigenvalue-ratio based on RMFA}
  \item{\code{"E_RM",}}{the rank-minimization based on IHR}
  \item{\code{"E_ER",}}{the eigenvalue-ratio based on IHR}
 }
}
  \item{max_iter}{
  Only if \code{method="E_RM"} or \code{method="E_ER"}, the maximum number of iterations in the iterative Huber regression algorithm. The default is 100.
}
  \item{c}{
  A constant to avoid vanishing denominators. The default is \eqn{10^{-4}}.
}
  \item{ep}{
  Only if \code{method="E_RM"} or \code{method="E_ER"}, the stopping critetion parameter in the iterative Huber regression algorithm. The default is \eqn{10^{-4} \times Tp_1 p_2}.
}
}
\details{
If \code{method="P"}, the number of factors \eqn{k_1} and \eqn{k_2} are estimated by \deqn{\hat{k}_1 = \arg \max_{j \leq k_{max}} \frac{\lambda _j (\bold{M}_c^w)}{\lambda _{j+1} (\bold{M}_c^w)}, \hat{k}_2 = \arg \max_{j \leq k_{max}} \frac{\lambda _j (\bold{M}_r^w)}{\lambda _{j+1} (\bold{M}_r^w)},} where \eqn{k_{max}} is a predetermined value larger than \eqn{k_1} and \eqn{k_2}. \eqn{\lambda _j(\cdot)} is the j-th largest eigenvalue of a nonnegative definitive matrix. See the function \code{\link{MHFA}} for the definition of \eqn{\bold{M}_c^w} and \eqn{\bold{M}_r^w}. For details, see He et al. (2023).

Define \eqn{D=\min({\sqrt{Tp_1}},\sqrt{Tp_2},\sqrt{p_1 p_2})},
\deqn{\hat{\bold{\Sigma}}_1=\frac{1}{T}\sum_{t=1}^T\hat{\bold{F}}_t \hat{\bold{F}}_t^\top, \hat{\bold{\Sigma}}_2=\frac{1}{T}\sum_{t=1}^T\hat{\bold{F}}_t^\top \hat{\bold{F}}_t,}
where \eqn{\hat{\bold{F}}_t, t=1, \dots, T} is estimated by IHR under the number of factor is \eqn{k_{max}}.

If \code{method="E_RM"}, the number of factors \eqn{k_1} and \eqn{k_2} are estimated by \deqn{\hat{k}_1=\sum_{i=1}^{k_{max}}I\left(\mathrm{diag}(\hat{\bold{\Sigma}}_1)>P_1\right), \hat{k}_2=\sum_{j=1}^{k_{max}}I\left(\mathrm{diag}(\hat{\bold{\Sigma}}_2) > P_2\right),}
where \eqn{I} is the indicator function. In practice, \eqn{P_1} is set as \eqn{\max \left(\mathrm{diag}(\hat{\bold{\Sigma}}_1)\right) \cdot D^{-2/3}}, \eqn{P_2} is set as \eqn{\max \left(\mathrm{diag}(\hat{\bold{\Sigma}}_2)\right) \cdot D^{-2/3}}.

If \code{method="E_ER"}, the number of factors \eqn{k_1} and \eqn{k_2} are estimated by \deqn{\hat{k}_1 = \arg \max_{i \leq k_{max}} \frac{\lambda _i (\hat{\bold{\Sigma}}_1)}{\lambda _{i+1} (\hat{\bold{\Sigma}}_1)+cD^{-2}}, \hat{k}_2 = \arg \max_{j \leq k_{max}} \frac{\lambda _j (\hat{\bold{\Sigma}}_2)}{\lambda _{j+1} (\hat{\bold{\Sigma}}_2)+cD^{-2}}.}
}
\value{
\item{\eqn{k_1}}{The estimated row factor number.}
\item{\eqn{k_2}}{The estimated column factor number.}
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
   
   KMHFA(X, kmax=6, method="P")
   \donttest{
   KMHFA(X, W1 = NULL, W2 = NULL, 6, "E_RM")
   KMHFA(X, W1 = NULL, W2 = NULL, 6, "E_ER")
   }
}