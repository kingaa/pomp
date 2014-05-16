\name{pomp-package}
\docType{package}
\alias{pomp-package}
\title{Partially-observed Markov processes}
\description{
  The \pkg{pomp} package provides facilities for inference on time series data using partially-observed Markov processes (AKA state-space models or nonlinear stochastic dynamical systems).
  One can use \pkg{pomp} to fit nonlinear, non-Gaussian dynamic models to time-series data.
  The first step in using \pkg{pomp} is to encode one's model and data in an object of class \code{pomp}.
  One does this via a call to \code{\link{pomp}}, which involves specifying the process and measurement components of the model in one or more of a variety of ways.
  Details on this are given in the documentation for the \code{\link{pomp}} function.
  Examples are given in the \sQuote{Introduction to \pkg{pomp}} document, in the demos (\code{demo(package='pomp')}), and via the \code{\link{pompExample}} function.

  \pkg{pomp} version \Sexpr[echo=F,stage=install,results=text]{packageDescription("pomp",fields="Version")} provides algorithms for
  \enumerate{
    \item simulation of stochastic dynamical systems; see \code{\link[=simulate-pomp]{simulate}}
    \item particle filtering (AKA sequential Monte Carlo or sequential importance sampling); see \code{\link{pfilter}}
    \item the iterated filtering method of Ionides et al. (2006, 2011); see \code{\link{mif}}
    \item the nonlinear forecasting algorithm of Kendall et al. (2005); see \code{\link{nlf}}
    \item the particle MCMC approach of Andrieu et al. (2010); see \code{\link{pmcmc}}
    \item the probe-matching method of Kendall et al. (1999, 2005); see \code{\link{probe.match}}
    \item a spectral probe-matching method (Reuman et al. 2006, 2008); see \code{\link{spect.match}}
    \item synthetic likelihood a la Wood (2010); see \code{\link{probe}}
    \item approximate Bayesian computation (Toni et al. 2009); see \code{\link{abc}}
    \item the approximate Bayesian sequential Monte Carlo scheme of Liu & West (2001); see \code{\link{bsmc}}
    \item basic trajectory matching; see \code{\link{traj.match}}.
  }
  The package also provides various tools for plotting and extracting information on models and data as well as an API for algorithm development.
  Future support for additional algorithms in envisioned.
  Much of the groundwork for \pkg{pomp} was laid by a working group of the National Center for Ecological Analysis and Synthesis (\acronym{NCEAS}), \dQuote{Inference for Mechanistic Models}.

  \pkg{pomp} is provided under the \acronym{GNU} Public License (\acronym{GPL}).
  Contributions are welcome, as are comments, feature requests, and bug reports.
  See the package website \url{http://pomp.r-forge.r-project.org} for more information, access to the package mailing list, links to the authors' websites, references to the literature, and up-to-date versions of the package source and documentation.
}
\section{Classes}{
  \pkg{pomp} makes extensive use of S4 classes.
  The basic class, \code{\link{pomp}}, encodes a partially-observed Markov process together with a uni- or multi-variate data set and (possibly) parameters.
}
\section{Documentation}{
  The \href{http://pomp.r-forge.r-project.org/vignettes/index.html}{\sQuote{Introduction to \pkg{pomp}}} document illustrates the facilities of the package using familiar stochastic processes.
  The \href{http://pomp.r-forge.r-project.org/vignettes/index.html}{\sQuote{Advanced topics in \pkg{pomp}}} document discusses the low-level interface and some more advanced modeling techniques.
  These documents can be viewed at \url{http://pomp.r-forge.r-project.org}.
}
\references{
  See the package website, \url{http://pomp.r-forge.r-project.org}, for the references.
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
