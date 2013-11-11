\name{logmeanexp}
\alias{logmeanexp}
\title{The log-mean-exp trick}
\description{
  \code{logmeanexp} computes the log-mean-exp of a set of numbers.
}
\usage{
logmeanexp(x)
}
\arguments{
  \item{x}{numeric}
}
\value{
  \code{log(mean(exp(x)))} computed so as to avoid over- or underflow.
}
\author{Aaron A. King \email{kingaa at umich dot edu}}
\examples{
  ## generate a bifurcation diagram for the Ricker map
  pompExample(ricker)
  ll <- replicate(n=5,logLik(pfilter(ricker,Np=1000)))
  ## an estimate of the log likelihood:
  logmeanexp(ll)
}
