\name{NLF}

\alias{NLF}

\title{Fit Model to Data Using Nonlinear Forecasting (NLF)}

\description{
 Calls an optimizer to maximize the nonlinear forecasting (NLF) goodness of fit, by
 simulating data from a model, fitting a nonlinear autoregressive
 model to the simulated time series (which may be multivariate) and
 using it to predict the data time series.
}

\usage{
 NLF(object, xstart, params, est, xfit, lags, nconverge=1000, nasymp=1000,
    seed = 1066, nrbf = 4, method = 'subplex', verbose = FALSE, gr = NULL, ...)
}

\arguments{
  \item{object}{A \code{pomp} object, with the data and model to fit to it.}
  \item{xstart}{Named numeric vector with guessed initial conditions.}
  \item{params}{Named numeric vector with guessed parameters.}
  \item{est}{Character vector containing the names of parameters to be
    estimated.}
  \item{xfit}{Character vector containing the names of state variables to be
    fit.}
 \item{lags}{A vector specifying the lags to use when constructing the
    nonlinear autoregressive prediction model. The first lag is the prediction
    interval.}
 \item{nconverge}{Number of convergence timesteps to be discarded
 from the model simulation}
 \item{nasymp}{Number of asymptotic timesteps to be recorded from the
 model simulation.}
 \item{seed}{Integer specifying the random number seed to use.  When fitting,
    it is usually best to always run the simulations with the same
    sequence of random numbers, which is accomplished by setting \code{seed} to an
    integer.  If you want a truly random simulation, set \code{seed=NULL}}
 \item{nrbf}{A scalar specifying the number of radial basis functions to be
    used at each lag.  See documentation to \code{\link{NLF.guts}} for details.}
 \item{method}{Optimization method.  Choices are \code{\link[subplex]{subplex}} and
 any of the methods used by \code{\link{optim}}}
 \item{verbose}{Logical; if \code{TRUE}, the negative log quasilikelihood and parameter
    values are printed at each iteration of the optimizer.}
   \item{gr}{Something for \code{optim}}
  \item{\dots}{Arguments that will be passed to \code{optim} in the
    \code{control} list.}

}

\details{

This is functionally a wrapper for \code{\link{NLF.objfun}}, which
does the statistical heavy lifting and should be consulted for
details.

}



\value{A list corresponding to the output from the optimizer, except
that the full parameter vector is returned (not just the ones
fitted), the LQL (and not -LQL) is reported, and xstart is
included.}

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
 Available online at
 \url{http://www2.bren.ucsb.edu/~kendall/pubs/1999Ecology.pdf}

 Kendall, B. E., Ellner, S. P., McCauley, E., Wood, S. N.,
 Briggs, C. J., Murdoch, W. W. and Turchin, P. (2005) Population cycles
 in the pine looper moth (\emph{Bupalus piniarius}): dynamical tests of
 mechanistic hypotheses. \emph{Ecological Monographs} \bold{75}, 259--276.
 Available online at
 \url{http://repositories.cdlib.org/postprints/818/}
}

\author{Stephen P. Ellner (spe2 at cornell dot edu) and Bruce E. Kendall
    (kendall at bren dot ucsb dot edu)}

\seealso{\code{\link{NLF.objfun}}, \code{\link{NLF.guts}}}

\examples{ }

\keyword{models}

\keyword{ts}
