\name{pomp-package}
\docType{package}
\alias{pomp-package}
\title{Partially-observed Markov processes}
\description{
  The \code{pomp} package provides facilities for inference on time series data using partially-observed Markov processes (AKA state-space models or nonlinear stochastic dynamical systems).
  The user encodes a model as a \code{\link{pomp}} object by providing functions specifying some or all of the model's process and measurement components.
  The package's algorithms for fitting models to data, simulating, etc. then call these functions.
  At the moment, algorithms are provided for
  particle filtering (AKA sequential Monte Carlo or sequential importance sampling, see \code{\link{pfilter}}),
  the likelihood maximization by iterated filtering (\acronym{MIF}) method of Ionides, Breto, and King (PNAS, 103:18438-18443, 2006, see \code{\link{mif}}),
  and the nonlinear forecasting algorithm of Kendall, Ellner, et al. (Ecol. Monog. 75:259-276, 2005, see \code{\link{nlf}}).
  Future support for a variety of other algorithms is envisioned.
  A working group of the National Center for Ecological Analysis and Synthesis (\acronym{NCEAS}), "Inference for Mechanistic Models", is currently implementing and testing additional methods for this package.

  The package is provided under the \acronym{GNU} Public License (\acronym{GPL}).
  Contributions are welcome, as are comments, suggestions for improvements, and bug reports.
}
\section{Classes}{
  The basic class, \code{\link{pomp}}, is provided to encode a partially-observed Markov process together with a uni- or multi-variate data set.
}
\section{Vignettes}{
  The vignette \sQuote{intro\_to\_pomp} illustrates the facilities of the package using familiar stochastic processes.
  Run \code{vignette("intro_to_pomp")} or look at the HTML documentation to view the vignette.
}
\author{Aaron A. King \email{kingaa at umich dot edu}}
\seealso{
  \code{\link{pomp}}, \code{\link{simulate-pomp}}, \code{\link{mif}}, \code{\link{nlf}}
}
\keyword{models}
\keyword{ts}
