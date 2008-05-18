\name{NLF.guts}

\alias{NLF.guts}

\title{Predict Data With Model Using Nonlinear Forecasting (NLF)}

\description{
 Uses a time series from a model simulation to construct a nonlinear
 auto-regressive model that predicts values at time \eqn{t} using values
 at past times.  This model is then applied to data time series, and
 the resulting sum of squared residuals is converted to a
 quasi-likelihood using assumptions of normality.
}

\usage{
 NLF.guts(data.mat, model.mat, lags, nrbf=4, verbose=FALSE, plotfit=FALSE,
    seas=FALSE, seas.data=NA, seas.sim=NA)
}

\arguments{
 \item{data.mat}{An array holding the data.  This array is of dimensions
    \code{nobs} x \code{ntdata}, where \code{nobs} is the number of
    observed variables and \code{ntdata} is the number of times at
    which observations were made.  Univariate data may be
    represented as a vector.  There can be no missing data, and the
    timesteps are assumed to be constant.}
 \item{model.mat}{An array holding the model simulation.  This array is of dimensions
    \code{nobs} x \code{ntmodel}, where \code{nobs} is the number of
    observed variables and \code{ntmodel} is the number of times at which
    observations were made.  Univariate data may be represented as a
    vector.  The number of observed variables in \code{model.mat}
    \emph{must} be the same as the number of observed variables in
    \code{data.mat}}
 \item{lags}{A vector specifying the lags to use when constructing the
    nonlinear autoregressive prediction model.  See Details for more
    information on how these should be constructed.}
 \item{nrbf}{A scalar specifying the number of radial basis functions to be
    used at each lag.  See Details.}
 \item{verbose}{Logical; if \code{TRUE}, various profiling and intermediate
    calculation values are printed.}
 \item{plotfit}{Logical; if \code{TRUE}, a diagnostic plot is produced.}
 \item{seas}{Logical; if \code{TRUE}, different models are fit for each
    season.  Currently not implemented, so has no effect.}
 \item{seas.data}{Not currently used (leave at default).}
 \item{seas.sim}{Not currently used (leave at default).}
}

\details{
 First consider the univariate case (scalar time series of data and
 model).  The following description is modified from Appendix B of
 Kendall et al. (2005).

The statistical forecasting model fitted to the model time series is
a generalized additive model of the form
 \deqn{
    x(t) = \sum_{j = 1}^n f\left[x(t - L_j) ;\theta_j \right] +
    \varepsilon_t
 }{
    x(t) = f[x(t-L_1); theta_1] + f[x(t-L_2); theta_2] + ... +
            f[x(t-L_n); theta_n] + epsilon(t)
 }
where \eqn{x(t)} is value of the state variable at time \eqn{t}
(where time is integer valued), the lags \eqn{(L_1, L_2, \dots,
L_n)}{(L_1, L_2, ..., L_n)} are integers specified in the calling
parameter \code{lags}, the ridge functions \eqn{f[x,\theta]}{f[x,
theta]} are each a sum of radial basis functions, and the process
noise variables \eqn{\varepsilon(t)}{epsilon(t)} are assumed to be
\emph{iid} Gaussian with zero mean.

The form of the ridge functions is
 \deqn{
    f[x;\theta] = \sum\limits_{k = 1}^r c_k \exp
    \left( - (x  - m_k )^2/2s^2 \right),
 }{
    f[x; theta] = c_1 exp\{-[(x - m_1)^2]/(2 s^2)\} + ... +
                   c_r exp\{-[(x - m_r)^2]/(2 s^2)\}
 }
so that the parameter vector for each ridge function is \eqn{\theta
= (s,c_1,m_1,c_2,m_2,\dots,c_r,m_r)}. That is, each \eqn{f} is the
sum of \eqn{r} (controlled by \code{nrbf} in the call) Gaussian
radial basis functions with centers at locations \eqn{m_k} and scale
parameter \eqn{s}. The values of \eqn{s} and of the \eqn{m_k} were
set based on the range \eqn{R} of the model time series. The
\eqn{m_k} were evenly spaced with the lowest value \eqn{R/10} below
the minimum of the the model time series and the largest value
\eqn{R/10} above the maximum of the the model time series.
\eqn{s=0.3R} \strong{(this may change)}. Given these, the full model
is linear in the remaining parameters (the \eqn{c}'s for each ridge
function) and therefore can be fitted very quickly by ordinary least
squares using \code{\link{lm}}.

The log quasi-likelihood, \eqn{L}, is defined as follows, based on
the assumption that the residuals (\eqn{\varepsilon}{epsilon}) are
normally distributed:
 \deqn{
    L = \sum_{t=1}^T -\frac{\varepsilon(t)^2}{\sigma^2} -
    \frac{1}{2} \log(2 \pi \sigma)
 }{
    L = [epsilon(1)^2 + ... + epsilon(T)^2]/(sigma^2) - log(2 pi
    sigma)/2.
 }
where \eqn{\sigma}{sigma} is the standard deviations of the
residuals and \eqn{T} is the number of fitted data points (which is
less than the total number of points in the data series by the
amount of the maximum lag).

When the data and model time series contain multiple state
variables, then each state variable is fit as a function of lagged
values of \emph{all} the state variables: if there are \eqn{q}
variables, then for state variable \eqn{i},
 \deqn{
    x_i(t) = \sum_{k=1}^q \sum_{j = 1}^n f\left[x_k(t - L_j) ;\theta_{jk} \right] +
    \varepsilon_{it}
 }{
    x_i(t) = f[x_1(t-L_1); theta_1,1] + f[x_1(t-L_2); theta_2,1] + ... +
            f[x_1(t-L_n); theta_n,1] + ... +
            f[x_q(t-L_1); theta_1,q] + f[x_q(t-L_2); theta_2,q] + ... +
            f[x_q(t-L_n); theta_n,q] +
            epsilon_i(t)
 }
The sums of squared residuals are converted to a quasi-likelihood
using the variance-covariance matrix of the residuals from each of
the state variables.

The background and rationale for NLF applied to univariate time
series can be found in Ellner et al (1998) and Kendall et al. (1999,
2005).}



\value{A scalar, the negative of the log quasi-likelihood.  Smaller
values indicate better predictions of the data.}

\references{

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

\note{This function will typically be called by
\code{\link{NLF.objfun}}; there is little error checking to make
sure that the dimensions of \code{data.mat} and \code{model.mat}
match appropriately, or that the values specified in \code{lags}
make sense.}

\author{Stephen P. Ellner (spe2 at cornell dot edu) and Bruce E. Kendall
    (kendall at bren dot ucsb dot edu)}

\seealso{\code{\link{NLF}}, \code{\link{NLF.objfun}}}

\examples{ }

\keyword{models}

\keyword{ts}
