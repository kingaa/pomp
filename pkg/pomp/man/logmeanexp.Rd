\name{logmeanexp}
\alias{logmeanexp}
\title{The log-mean-exp trick}
\description{
  \code{logmeanexp} computes the log-mean-exp of a set of numbers.
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
\author{Aaron A. King \email{kingaa at umich dot edu}}
\examples{
  ## generate a bifurcation diagram for the Ricker map
  pompExample(ricker)
  ll <- replicate(n=5,logLik(pfilter(ricker,Np=1000)))
  ## an estimate of the log likelihood:
  logmeanexp(ll)
  logmeanexp(ll,se=TRUE)
}
