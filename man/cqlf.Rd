\name{cqlf}
\alias{cqlf}
\title{Composite quasilikelihood function.}
\description{
Compute the quasilikelihood function given a set of simulated data and an actual data set of the same length.
}
\usage{
cqlf(simdata, data, maxlag=ncol(data))
}
\arguments{
  \item{simdata}{A matrix containing the simulated data sets, each data set in its own column.}
  \item{data}{The actual data set.  This must be a vector or 1-D array.}
  \item{maxlag}{The maximum lag to be used in the computation.}
}
\details{
  Suppose the data are \eqn{\{y_t\}_{t=1}^T} and that, for given values of the parameters, the model predicts mean and variance \eqn{\mu_t}, \eqn{\sigma_t^2} in \eqn{y_t}.
  Let \eqn{z_t=((y_t-\mu_t)/\sigma_t)} and \eqn{\rho_{st}} be the expected correlation between times \eqn{s} and \eqn{t}.
  The composite quasilikelihood \eqn{Q} is defined by \deqn{-2\log{Q}=\sum_{t=1}^{T}\sum_{s=t-m}^{t-1} \{(z_s^2+z_t^2-2 \rho_{st} z_s z_t)/(1-\rho_{st}^2)+\log{\sigma_s}+\log{\sigma_t}+\log{(1-\rho_{st}^2)}\}}.
}
\value{
\eqn{-2 \log{Q}}
}
\references{
  S. R. Lele (2006)
  Sampling variability and estimates of density dependence: a composite-likelihood approach.
  Ecology 87:189-202

  S. R. Lele and M. L. Taper (2002)
  A composite likelihood approach to (co)variance components estimation.
  Journal of Statistical Planning and Inference 103:117-135.

  F. Curriero and S. Lele (1999)
  A composite likelihood approach to semivariogram estimation.
  Journal of Agricultural, Biological and Environmental Statistics 4:9-28.

  P. Heagerty and S. R. Lele (1998) A composite likelihood approach to binary data in space.
  Journal of American Statistical Association 93:1099-1111.
}
\author{Aaron A. King (kingaa at umich dot edu)}
\examples{}
\keyword{models}
\keyword{ts}
