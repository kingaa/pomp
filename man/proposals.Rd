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
  mvn.rw(rw.var)
  mvn.diag.rw(rw.sd)
  mvn.rw.adaptive(rw.sd, rw.var, scale.start = NA, scale.cooling = 0.999,
                  shape.start = NA, target = 0.234, max.scaling = 50)
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
  \item{scale.start, scale.cooling, shape.start, target, max.scaling}{
    parameters to control the proposal adaptation algorithm.
    Beginning with MCMC iteration \code{scale.start}, the scale of the proposal covariance matrix will be adjusted in an effort to match the \code{target} acceptance ratio.
    This initial scale adjustment is "cooled", that is the adjustment diminishes as the chain moves along.
    The parameter \code{scale.cooling} specifies the cooling schedule: at n iterations after \code{scale.start}, the current scaling factor is multiplied with \code{scale.cooling}^n.
    A maximum scaling factor can be specified with \code{max.scaling}.
    After \code{shape.start} accepted proposal, a scaled empirical covariance matrix will be used for the proposals according to Roberts and Rosenthal (2009).
  }
}
\value{
  Each of these constructs a function suitable for use as the \code{proposal} argument of \code{pmcmc} or \code{abc}.
  Given a parameter vector, each such function returns a single draw from the corresponding proposal distribution.
}
\seealso{\code{\link{pmcmc}}, \code{\link{abc}}}
\references{
  Gareth O. Roberts and Jeffrey S. Rosenthal.
  Examples of Adaptive MCMC.
  J. Comput. Graph. Stat., 18:349--367, 2009.
}
\author{Aaron A. King, Sebastian Funk}
\keyword{programming}
