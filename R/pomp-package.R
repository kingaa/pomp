##' Inference for partially observed Markov processes
##'
##' The \pkg{pomp} package provides facilities for inference on time series
##' data using partially-observed Markov process (\acronym{POMP}) models.
##' These models are also known as state-space models, hidden Markov models, or
##' nonlinear stochastic dynamical systems.  One can use \pkg{pomp} to fit
##' nonlinear, non-Gaussian dynamic models to time-series data.  The package is
##' both a set of tools for data analysis and a platform upon which statistical
##' inference methods for \acronym{POMP} models can be implemented.
##'
##' @name pomp-package
##' @aliases pomp,package pomp-package
##' @docType package
##' @author Aaron A. King
##' @family information on model implementation
##' @family \pkg{pomp} parameter estimation methods
##' @family elementary POMP methods
##' @keywords ts models multivariate
##'
##' @section Data analysis using \pkg{pomp}:
##' \pkg{pomp} provides algorithms for
##' \enumerate{
##' \item simulation of stochastic
##' dynamical systems; see \code{\link[=simulate-pomp]{simulate}}
##' \item
##' particle filtering (AKA sequential Monte Carlo or sequential importance
##' sampling); see \code{\link{pfilter}}
##' \item the iterated filtering methods
##' of Ionides et al. (2006, 2011, 2015); see \code{\link{mif2}}
##' \item the
##' nonlinear forecasting algorithm of Kendall et al. (2005); see
##' \code{\link{nlf}}
##' \item the particle MCMC approach of Andrieu et al. (2010); see \code{\link{pmcmc}}
##' \item the probe-matching method of Kendall et al. (1999, 2005); see \code{\link{probe.match}}
##' \item a spectral probe-matching method (Reuman et al. 2006, 2008); see
##' \code{\link{spect.match}}
##' \item synthetic likelihood a la Wood (2010); see \code{\link{probe}}
##' \item approximate Bayesian computation (Toni et al. 2009); see \code{\link{abc}}
##' \item the approximate Bayesian sequential
##' Monte Carlo scheme of Liu & West (2001); see \code{\link{bsmc2}}
##' \item ensemble and ensemble adjusted Kalman filters; see \code{\link{kalman}}
##' \item simple trajectory matching; see \code{\link{traj.match}}.
##' }
##' The package
##' also provides various tools for plotting and extracting information on
##' models and data.
##'
##' @references
##' A. A. King, D. Nguyen, and E. L. Ionides (2016) Statistical
##' Inference for Partially Observed Markov Processes via the Package
##' \pkg{pomp}.  \emph{Journal of Statistical Software} 69(12): 1--43.  An
##' updated version of this paper is available on the
##' \href{https://kingaa.github.io/pomp/docs.html}{package website}.
##'
##' See the package website, \url{https://kingaa.github.io/pomp/}, for more
##' references.
##'
##' @useDynLib pomp, .registration = TRUE, .fixes="P_"
##' @import methods
NULL
