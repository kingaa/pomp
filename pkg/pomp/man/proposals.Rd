\name{MCMC proposal distributions}
\alias{MCMC proposal distributions}
\alias{MCMC proposal functions}
\alias{mvn.rw}
\alias{mvn.diag.rw}
\title{MCMC proposal distributions}
\description{
  Functions to construct proposal distributions for use with MCMC methods.
}
\usage{
  mvn.diag.rw(rw.sd)
  mvn.rw(rw.var)
}
\arguments{
  \item{rw.sd}{
    named numeric vector;
    random-walk SDs for a multivariate normal random-walk proposal with diagonal variance-covariance matrix.
  }
  \item{rw.var}{
    square numeric matrix with row- and column-names.
    Specifies the variance-covariance matrix for a multivariate normal random-walk proposal distribution.
  }
}
\value{
  Functions taking a single argument.
  Given a parameter vector, each returns a single draw from the corresponding proposal distribution.
}
\seealso{\code{\link{pmcmc}}, \code{\link{abc}}}
\author{Aaron A. King \email{kingaa at umich dot edu}}
\keyword{models}
