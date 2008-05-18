\name{NLF.objfun}

\alias{NLF.objfun}

\title{Evaluate Fit of Model to Data Using Nonlinear Forecasting (NLF)}

\description{
 Evaluates the nonlinear forecasting (NLF) goodness of fit, by
 simulating data from a model, fitting a nonlinear autoregressive
 model to the simulated time series (which may be multivariate) and
 using it to predict the data time series.
}

\usage{
 NLF.objfun(params.fitted, object, params, par.index, ts.index, xstart,
    times, lags, seed, nrbf=4, verbose=FALSE)
}

\arguments{
 \item{params.fitted}{Numeric vector (of length \code{nfit})
    with the subset of the parameters to be
    fitted (these are kept separately so that they can be adjusted by an
    optimizer)}
 \item{object}{A \code{pomp} object, with the data and model to fit to it}
 \item{params}{Numeric vector with parameter values}
 \item{par.index}{Integer vector indicating of length \code{nfit} indicating
    which elements of \code{params} should be replaced by the values in
    \code{params.fitted}}
 \item{ts.index}{Integer vector specifying which rows (state variables) of
    the \code{data} element of \code{object} should be fitted}
 \item{xstart}{Numeric vector of initial conditions}
 \item{times}{Numeric vector of times at which to evaluate the model.  It is
    recommended to drop a largish number of convergence iterations; do this
    by leaving a large gap between the first (initial time) and second
    (first reported time) value of \code{times}}
 \item{lags}{A vector specifying the lags to use when constructing the
    nonlinear autoregressive prediction model}
 \item{seed}{Integer specifying the random number seed to use.  When fitting,
    it is usually best to always run the simulations with the same
    sequence of random numbers, which is accomplished by setting \code{seed} to an
    integer.  If you want a truly random simulation, set \code{seed=NULL}}
 \item{nrbf}{A scalar specifying the number of radial basis functions to be
    used at each lag}
 \item{verbose}{Logical; if \code{TRUE}, the return value and parameter
    values are printed.  This is useful to track progress when it is being
    called by an optimizer.}
}

\details{

This is functionally a wrapper for \code{\link{NLF.guts}}, which
does the statistical heavy lifting and should be consulted for
details.  It is structured to be called from an optimizer (such as
\code{\link[subplex]{subplex}} or \code{\link{optim}}.

}



\value{A scalar, the negative of the log quasi-likelihood.  Smaller
values indicate better predictions of the data.}

\references{
 The following papers describe and motivate the NLF approach to model
 fitting:

 Ellner, S. P., Bailey, B. A., Bobashev, G. V., Gallant, A. R.,
 Grenfell, B. T. and Nychka D. W. (1998)  Noise and nonlinearity in
 measles epidemics: combining mechanistic and statistical approaches
 to population modeling. \emph{American Naturalist} \bold{151}, 425--440.

 Kendall, B. E., Briggs, C. J., Murdoch, W. W., Turchin, P.,
 Ellner, S. P., McCauley, E., Nisbet, R. M. and Wood S. N. (1999)
 Why do populations cycle? A synthesis of statistical and mechanistic
 modeling approaches.  \emph{Ecology} \bold{80}, 1789--1805.

 Kendall, B. E., Ellner, S. P., McCauley, E., Wood, S. N.,
 Briggs, C. J., Murdoch, W. W. and Turchin, P. (2005) Population cycles
 in the pine looper moth (\emph{Bupalus piniarius}): dynamical tests of
 mechanistic hypotheses. \emph{Ecological Monographs} \bold{75}, 259--276.
}

\author{Stephen P. Ellner (spe2 at cornell dot edu) and Bruce E. Kendall
    (kendall at bren dot ucsb dot edu)}

\seealso{\code{\link{NLF}}, \code{\link{NLF.guts}}}

\examples{ }

\keyword{models}

\keyword{ts}
