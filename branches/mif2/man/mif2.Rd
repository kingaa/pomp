\name{mif2}
\docType{methods}
\alias{mif2}
\alias{mif2,mif2d.pomp-method}
\alias{mif2-mif2d.pomp}
\alias{mif2,pfilterd.pomp-method}
\alias{mif2-pfilterd.pomp}
\alias{mif2,pomp-method}
\alias{mif2-pomp}
\alias{continue,mif2d.pomp-method}
\alias{continue-mif2}
\title{The MIF2 algorithm}
\description{
  The MIF2 algorithm for estimating the parameters of a partially-observed Markov process.
}
\usage{
\S4method{mif2}{pomp}(object, Nmif = 1, start, Np, perturb.fn,
    tol = 1e-17, max.fail = Inf, verbose = getOption("verbose"),
    transform = FALSE, \dots)
\S4method{mif2}{pfilterd.pomp}(object, Nmif = 1, Np, tol, \dots)
\S4method{mif2}{mif2d.pomp}(object, Nmif, start, Np, perturb.fn,
    tol, transform, \dots)
\S4method{continue}{mif2d.pomp}(object, Nmif = 1, \dots)
}
\arguments{
  \item{object}{
    An object of class \code{pomp}.
  }
  \item{Nmif}{
    The number of MIF iterations to perform.
  }
  \item{start}{
    named numerical vector;
    the starting guess of the parameters.
  }
  \item{Np}{
    the number of particles to use in filtering.
    This may be specified as a single positive integer, in which case the same number of particles will be used at each timestep.
    Alternatively, if one wishes the number of particles to vary across timestep, one may specify \code{Np} either as a vector of positive integers (of length \code{length(time(object,t0=TRUE))}) or as a function taking a positive integer argument.
    In the latter case, \code{Np(k)} must be a single positive integer, representing the number of particles to be used at the \code{k}-th timestep:
    \code{Np(0)} is the number of particles to use going from \code{timezero(object)} to \code{time(object)[1]},
    \code{Np(1)}, from \code{timezero(object)} to \code{time(object)[1]},
    and so on, while when \code{T=length(time(object,t0=TRUE))},
    \code{Np(T)} is the number of particles to sample at the end of the time-series.
  }
  \item{tol}{
    See the description under \code{\link{pfilter}}.
  }
  \item{max.fail}{
    See the description under \code{\link{pfilter}}.
  }
  \item{verbose}{
    logical; if TRUE, print progress reports.
  }
  \item{transform}{
    logical;
    if \code{TRUE}, optimization is performed on the transformed scale.
  }
  \item{\dots}{
    additional arguments that override the defaults.
  }
}
\section{Re-running MIF Iterations}{
}
\section{Continuing MIF Iterations}{
}
\section{Using MIF to estimate initial-value parameters only}{
}
\section{Details}{
}
\references{
  E. L. Ionides, A. Bhadra, Y. Atchad{\\'e}, & A. A. King,
  Iterated filtering,
  Annals of Statistics, 39:1776--1802, 2011.
}
\author{Aaron A. King \email{kingaa at umich dot edu}}
\seealso{
  \code{\link{pomp}}, \code{\link{pomp-class}}, \code{\link{pfilter}}.
}
\keyword{ts}
