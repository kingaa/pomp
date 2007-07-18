\name{pomp-package}
\docType{package}
\alias{pomp-package}

\title{Partially-observed Markov processes}
\description{The `pomp' package provides facilities for
  partially-observed Markov processes.}
\section{Classes}{
  The basic class, \code{pomp}, is provided to encode a
  partially-observed Markov process together with a multivariate data
  set.

  The class \code{mif} derives from class \code{pomp} and encodes the
  results of fitting the model to the data by the MIF algorithm.
}
\section{Vignettes}{
  The vignette `random\_walk\_example' illustrates the facilities of the
  package using simple random-walk processes.  Run
  \code{vignette('random_walk_example')} or look at the HTML
  documentation to view the vignette.
}
\references{}
\author{Aaron A. King (kingaa at umich dot edu)}
\seealso{\code{\link{pomp-class}}, \code{\link{pomp}},
  \code{\link{mif-class}}, \code{\link{mif}}
}
\keyword{models}
\keyword{ts}
