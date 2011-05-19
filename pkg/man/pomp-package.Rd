\name{pomp-package}
\docType{package}
\alias{pomp-package}
\title{Partially-observed Markov processes}
\description{
%%  This is documentation for version \Sexpr[echo=F,stage=install,results=text]{packageDescription("pomp",fields="Version")} of \pkg{pomp}.

  The \pkg{pomp} package provides facilities for inference on time series data using partially-observed Markov processes (AKA state-space models or nonlinear stochastic dynamical systems).
  One might use \pkg{pomp} to fit a model to time-series data or simply to simulate a stochastic dynamical model.
  The first step in using \pkg{pomp} is to encode one's model and data in an object of class \code{pomp}.
  One does this via a call to \code{\link{pomp}}, which involves specifying the process and measurement components of the model in one or more of a variety of ways.
  Details on this are given in the documentation for the \code{\link{pomp}} function and examples are given in the \sQuote{intro_to_pomp} vignette.

  Currently, \pkg{pomp} provides algorithms for
  (i) simulation of stochastic dynamical systems (see \code{\link[=simulate-pomp]{simulate}}),
  (ii) particle filtering (AKA sequential Monte Carlo or sequential importance sampling, see \code{\link{pfilter}}),
  (iii) the likelihood maximization by iterated filtering (\acronym{MIF}) method of Ionides, Breto, and King (PNAS, 103:18438-18443, 2006, see \code{\link{mif}}),
  and (iv) the nonlinear forecasting algorithm of Kendall, Ellner, et al. (Ecol. Monog. 75:259-276, 2005, see \code{\link{nlf}}).
  Future support for additional algorithms in envisioned, and implementations of
  the Bayesian sequential Monte Carlo approach of Liu \& West and generalized probe-matching methods are currently under development.
  Much of the work in \pkg{pomp} has been done under the auspices of a working group of the National Center for Ecological Analysis and Synthesis (\acronym{NCEAS}), "Inference for Mechanistic Models".

  The package is provided under the \acronym{GNU} Public License (\acronym{GPL}).
  Contributions are welcome, as are comments, suggestions for improvements, and bug reports.
}
\section{Classes}{
  \pkg{pomp} makes extensive use of S4 classes.
  The basic class, \code{\link{pomp}}, encodes a partially-observed Markov process together with a uni- or multi-variate data set and (possibly) parameters.
}
\section{Vignettes}{
  The vignette \sQuote{intro_to_pomp} illustrates the facilities of the package using familiar stochastic processes.
  Run \code{vignette("intro_to_pomp")} or look at the HTML documentation to view the vignette.
}
\author{Aaron A. King \email{kingaa at umich dot edu}}
\seealso{
  \code{\link{pomp}},
  \code{\link{pfilter}},
  \code{\link[=simulate-pomp]{simulate}},
  \code{\link{trajectory}},
  \code{\link{mif}},
  \code{\link{nlf}},
  \code{\link{probe.match}},
  \code{\link{traj.match}},
  \code{\link{bsmc}},
  \code{\link{pmcmc}}
}
\keyword{models}
\keyword{ts}
