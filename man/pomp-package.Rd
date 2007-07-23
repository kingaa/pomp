\name{pomp-package}
\docType{package}
\alias{pomp-package}
\title{Partially-observed Markov processes}
\description{The `pomp' package provides facilities for inference using
  partially-observed Markov processes  (AKA state-space models or
  nonlinear stochastic dynamical systems). The user provides functions
  specifying the model's process and measurement components. The
  package's algorithms are built on top of these functions. At the
  moment, algorithms are provided for particle filtering (AKA sequential
  Monte Carlo or sequential importance sampling) and the likelihood
  maximization by iterated filtering (MIF) method of Ionides, Breto, and
  King (PNAS, 103:18438-18443, 2006). Future support for a variety of
  other algorithms is envisioned. A working group of the National Center
  for Ecological Analysis and Synthesis (NCEAS), "Inference for
  Mechanistic Models", is currently implementing additional methods for
  this package.

  The package is provided under the GPL. Contributions are welcome, as
  are comments, suggestions, and bug reports.
}
\section{Classes}{
  The basic class, \code{pomp}, is provided to encode a
  partially-observed Markov process together with a multivariate data
  set.

  The class \code{mif} derives from class \code{pomp} and encodes the
  results of fitting the model to the data by the MIF algorithm.
}
\section{Vignettes}{
  The vignette `random\_walk' illustrates the facilities of the
  package using simple random-walk processes.  Run
  \code{vignette('random_walk')} or look at the HTML
  documentation to view the vignette.
}
\author{Aaron A. King (kingaa at umich dot edu)}
\seealso{\code{\link{pomp-class}}, \code{\link{pomp}},
  \code{\link{mif-class}}, \code{\link{mif}}
}
\keyword{models}
\keyword{ts}
