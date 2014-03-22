\name{pomp-package}
\docType{package}
\alias{pomp-package}
\title{Partially-observed Markov processes}
\description{
%%  This is documentation for version \Sexpr[echo=F,stage=install,results=text]{packageDescription("pomp",fields="Version")} of \pkg{pomp}.

  The \pkg{pomp} package provides facilities for inference on time series data using partially-observed Markov processes (AKA state-space models or nonlinear stochastic dynamical systems).
  One can use \pkg{pomp} to fit nonlinear, non-Gaussian dynamic models to time-series data.
  The first step in using \pkg{pomp} is to encode one's model and data in an object of class \code{pomp}.
  One does this via a call to \code{\link{pomp}}, which involves specifying the process and measurement components of the model in one or more of a variety of ways.
  Details on this are given in the documentation for the \code{\link{pomp}} function and examples are given in the \sQuote{intro_to_pomp} vignette.

  Currently, \pkg{pomp} provides algorithms for
  (i) simulation of stochastic dynamical systems (see \code{\link[=simulate-pomp]{simulate}}),
  (ii) particle filtering (AKA sequential Monte Carlo or sequential importance sampling), see \code{\link{pfilter}}),
  (iii) the iterated filtering method of Ionides et al. (2006), see \code{\link{mif}}),
  (iv) the nonlinear forecasting algorithm of Kendall et al. (2005), see \code{\link{nlf}}),
  (v) the particle MCMC approach of Andrieu et al. (2010), see \code{\link{pmcmc}},
  (vi) basic trajectory matching, see \code{\link{traj.match}},
  (vi) the probe-matching method of Wood (2010) and Kendall et al. (1999), see \code{\link{probe.match}},
  (vii) a spectral probe-matching method (Reuman et al., 2006), see \code{\link{spect.match}}.
  See the package website \url{http://pomp.r-forge.r-project.org} for these references.
  The package also provides various tools for plotting and extracting information on models and data as well as an API for algorithm development.
  Future support for additional algorithms in envisioned, and implementations of
  the Bayesian sequential Monte Carlo approach of Liu \& West.
  Much of the work in \pkg{pomp} has been done under the auspices of a working group of the National Center for Ecological Analysis and Synthesis (\acronym{NCEAS}), "Inference for Mechanistic Models".

  The package is provided under the \acronym{GNU} Public License (\acronym{GPL}).
  Contributions are welcome, as are comments, suggestions for improvements, and bug reports.
  See the package website \url{http://pomp.r-forge.r-project.org} for more information, access to the package mailing list, links to the authors' websites, and references to the literature.
}
\section{Classes}{
  \pkg{pomp} makes extensive use of S4 classes.
  The basic class, \code{\link{pomp}}, encodes a partially-observed Markov process together with a uni- or multi-variate data set and (possibly) parameters.
}
\section{Vignettes}{
  The vignette \sQuote{Introduction to pomp} illustrates the facilities of the package using familiar stochastic processes.
  Methods for accelerating your codes are discussed in the \sQuote{Advanced topics in pomp} vignette.
  These vignettes are installed in the package \code{doc} directory and are accessible via the HTML documentation.
  They can also be viewed at \url{http://pomp.r-forge.r-project.org}.
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
