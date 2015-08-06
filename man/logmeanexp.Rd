\name{logmeanexp}
\alias{logmeanexp}
\title{The log-mean-exp trick}
\description{
  \code{logmeanexp} computes \deqn{\log\frac{1}{N}\sum_{n=1}^N\!e^x_i,}{log mean exp(x_i),} avoiding over- and under-flow in doing so.
  It can optionally return an estimate of the standard error in this quantity.
}
\usage{
logmeanexp(x, se = FALSE)
}
\arguments{
  \item{x}{numeric}
  \item{se}{logical; give approximate standard error?}
}
\value{
  \code{log(mean(exp(x)))} computed so as to avoid over- or underflow.
  If \code{se = FALSE}, the approximate standard error is returned as well.
}
\details{
  When \code{se = TRUE}, \code{logmeanexp} uses a jackknife estimate of the variance in \eqn{log(x)}.
}
\author{Aaron A. King}
\examples{
  ## generate a bifurcation diagram for the Ricker map
  pompExample(ricker)
  ll <- replicate(n=5,logLik(pfilter(ricker,Np=1000)))
  ## an estimate of the log likelihood:
  logmeanexp(ll)
  ## with standard error:
  logmeanexp(ll,se=TRUE)
}
