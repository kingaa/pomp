\name{MCMC proposal distributions}
\alias{MCMC proposal distributions}
\alias{MCMC proposal functions}
\alias{mvn.rw}
\alias{mvn.diag.rw}
\alias{mvn.rw.adaptive}
\title{MCMC proposal distributions}
\description{
  Functions to construct proposal distributions for use with MCMC methods.
}
\usage{
  mvn.rw(rw.var, log = FALSE)
  mvn.diag.rw(rw.sd, log = FALSE)
  mvn.rw.adaptive(rw.sd, rw.var, scale.start = NA, scale.cooling = 0.999,
                  shape.start = NA, target = 0.234, max.scaling = 50, \dots)
}
\arguments{
  \item{rw.var}{
    square numeric matrix with row- and column-names.
    Specifies the variance-covariance matrix for a multivariate normal random-walk proposal distribution.
  }
  \item{rw.sd}{
    named numeric vector;
    random-walk SDs for a multivariate normal random-walk proposal with diagonal variance-covariance matrix.
  }
  \item{log}{
    logical; if \code{TRUE}, a log-normal proposal distribution is used.
  }
  \item{scale.start, scale.cooling, shape.start, target, max.scaling}{
    parameters to control the proposal adaptation algorithm.
    Beginning with MCMC iteration \code{scale.start}, the scale of the proposal covariance matrix will be adjusted in an effort to match the \code{target} acceptance ratio.
    The parameters \code{scale.cooling} and \code{max.scaling} adjust the scale adaptation.
    Beginning with iteration \code{shape.start}, a scaled empirical covariance matrix will be used for the proposals.
  }
}
\value{
  Functions taking a single argument.
  Given a parameter vector, each returns a single draw from the corresponding proposal distribution.
}
\seealso{\code{\link{pmcmc}}, \code{\link{abc}}}
\author{Aaron A. King \email{kingaa at umich dot edu}}
\keyword{programming}
