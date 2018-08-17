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
##'
##' @name pomp-package
##' @docType package
##'
##' @section Data analysis using \pkg{pomp}:
##'
##' \pkg{pomp} version \Sexpr{packageDescription("pomp",fields="Version")}
##' provides algorithms for
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
##' \item ensemble and ensemble adjusted Kalman filters; see \code{\link{enkf}}
##' \item simple trajectory matching; see \code{\link{traj.match}}.
##' }
##' The package
##' also provides various tools for plotting and extracting information on
##' models and data.
##'
##' @author Aaron A. King
##'
##' @seealso \code{\link{pfilter}}, \code{\link[=simulate-pomp]{simulate}},
##' \code{\link{mif2}}, \code{\link{nlf}}, \code{\link{probe}},
##' \code{\link{traj.match}}, \code{\link{bsmc2}}, \code{\link{pmcmc}}
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
NULL

##' Model for Nicholson's blowflies.
##'
##' \code{blowflies1} and \code{blowflies2} are \sQuote{pomp} objects encoding
##' stochastic delay-difference models.
##'
##' The data are from "population I", a control culture in one of A. J.
##' Nicholson's experiments with the Australian sheep-blowfly \emph{Lucilia
##' cuprina}.  The experiment is described on pp. 163--4 of Nicholson (1957).
##' Unlimited quantities of larval food were provided; the adult food supply
##' (ground liver) was constant at 0.4g per day.  The data were taken from the
##' table provided by Brillinger et al. (1980).
##'
##' The models are discrete delay equations: \deqn{R(t+1) \sim
##' \mathrm{Poisson}(P N(t-\tau) \exp{(-N(t-\tau)/N_{0})} e(t+1) {\Delta}t)}{%
##' R[t+1] ~ Poisson(P N[t-tau] exp(-N[t-tau]/N0) e[t+1] dt)} \deqn{S(t+1) \sim
##' \mathrm{binomial}(N(t),\exp{(-\delta \epsilon(t+1) {\Delta}t)})}{% S[t+1] ~
##' binomial(N[t],exp(-delta eps[t+1] dt))} \deqn{N(t) =
##' R(t)+S(t)}{N[t]=R[t]+S[t]} where \eqn{e(t)}{e[t]} and
##' \eqn{\epsilon(t)}{eps[t]} are Gamma-distributed i.i.d. random variables
##' with mean 1 and variances \eqn{{\sigma_p^2}/{{\Delta}t}}{sigma.p^2/dt},
##' \eqn{{\sigma_d^2}/{{\Delta}t}}{sigma.d^2/dt}, respectively.
##' \code{blowflies1} has a timestep (\eqn{{\Delta}t}{dt}) of 1 day, and
##' \code{blowflies2} has a timestep of 2 days.  The process model in
##' \code{blowflies1} thus corresponds exactly to that studied by Wood (2010).
##' The measurement model in both cases is taken to be \deqn{y(t) \sim
##' \mathrm{negbin}(N(t),1/\sigma_y^2)}{y[t] ~ negbin(N[t],1/sigma.y^2)}, i.e.,
##' the observations are assumed to be negative-binomially distributed with
##' mean \eqn{N(t)}{N[t]} and variance \eqn{N(t)+(\sigma_y
##' N(t))^2}{N[t]+(sigma.y N[t])^2}.
##'
##' Do \preformatted{pompExample(blowflies,show=TRUE)} to view the code that
##' constructs these pomp objects.
##'
##' @name blowflies
##' @aliases blowflies blowflies1 blowflies2
##' @docType data
##' @seealso \code{\link{pomp}}
##' @references
##' A. J. Nicholson (1957)
##' The self-adjustment of populations to change.
##' Cold Spring Harbor Symposia on Quantitative Biology,
##' \bold{22}, 153--173.
##'
##' Y. Xia and H. Tong (2011) Feature Matching in Time Series Modeling.
##' \emph{Statistical Science} \bold{26}, 21--46.
##'
##' E. L. Ionides (2011) Discussion of ``Feature Matching in Time Series
##' Modeling'' by Y. Xia and H. Tong.  \emph{Statistical Science} \bold{26},
##' 49--52.
##'
##' S. N. Wood (2010) Statistical inference for noisy nonlinear ecological
##' dynamic systems.  \emph{Nature} \bold{466}, 1102--1104.
##'
##' W. S. C. Gurney, S. P. Blythe, and R. M. Nisbet (1980) Nicholson's
##' blowflies revisited.  \emph{Nature} \bold{287}, 17--21.
##'
##' D. R. Brillinger, J. Guckenheimer, P. Guttorp and G. Oster (1980) Empirical
##' modelling of population time series: The case of age and density dependent
##' rates.  in G. Oster (ed.), Some Questions in Mathematical Biology, vol. 13,
##' pp. 65--90.  American Mathematical Society, Providence.
##' @keywords models datasets
##' @examples
##'
##' pompExample(blowflies)
##' plot(blowflies1)
##' plot(blowflies2)
##'
NULL

##' Model of cholera transmission for historic Bengal.
##'
##' \code{dacca} is a \sQuote{pomp} object containing census and cholera
##' mortality data from the Dacca district of the former British province of
##' Bengal over the years 1891 to 1940 together with a stochastic differential
##' equation transmission model.  The model is that of King et al. (2008).  The
##' parameters are the MLE for the SIRS model with seasonal reservoir.
##'
##' Data are provided courtesy of Dr. Menno J. Bouma, London School of Tropical
##' Medicine and Hygiene.
##'
##' \code{dacca} is a \sQuote{pomp} object containing the model, data, and MLE
##' parameters.  Parameters that naturally range over the positive reals are
##' log-transformed; parameters that range over the unit interval are
##' logit-transformed; parameters that are naturally unbounded or take integer
##' values are not transformed.
##'
##' @name dacca
##' @docType data
##' @seealso \code{\link{sir}}, \code{\link{pomp}}
##' @references A. A. King, E. L. Ionides, M. Pascual, and M. J. Bouma,
##' Inapparent infections and cholera dynamics, Nature, 454:877-880, 2008
##' @keywords models datasets
##' @examples
##'
##' pompExample(dacca)
##' plot(dacca)
##' #MLEs on the natural scale
##' coef(dacca)
##' #MLEs on the transformed scale
##' coef(dacca,transform=TRUE)
##' plot(simulate(dacca))
##'
NULL

##' Gompertz model with log-normal observations.
##'
##' \code{gompertz} is a \sQuote{pomp} object encoding a stochastic Gompertz
##' population model with log-normal measurement error.
##'
##' The state process is \eqn{X_{t+1} = K^{1-S} X_{t}^S
##' \epsilon_{t}}{X[t+1]=K^(1-S) X[t]^S eps[t]}, where \eqn{S=e^{-r}}{S=e^{-r}}
##' and the \eqn{\epsilon_t}{eps[t]} are i.i.d. lognormal random deviates with
##' variance \eqn{\sigma^2}{sigma^2}.  The observed variables \eqn{Y_t} are
##' distributed as
##' \eqn{\mathrm{lognormal}(\log{X_t},\tau)}{lognormal(log(X[t]),tau)}.
##' Parameters include the per-capita growth rate \eqn{r}, the carrying
##' capacity \eqn{K}, the process noise s.d. \eqn{\sigma}{sigma}, the
##' measurement error s.d. \eqn{\tau}{tau}, and the initial condition
##' \eqn{X_0}{X[0]}.  The \sQuote{pomp} object includes parameter
##' transformations that log-transform the parameters for estimation purposes.
##'
##' @name gompertz
##' @docType data
##' @seealso \code{pomp}, \code{ricker}, and the tutorials at
##' \url{https://kingaa.github.io/pomp/}.
##' @keywords models datasets
##' @examples
##'
##' pompExample(gompertz)
##' plot(gompertz)
##' coef(gompertz)
##' coef(gompertz,transform=TRUE)
##'
NULL

##' Historical childhood disease incidence data
##'
##' \code{LondonYorke} is a data frame containing the monthly number of
##' reported cases of chickenpox, measles, and mumps from two American cities
##' (Baltimore and New York) in the mid-20th century (1928--1972).
##'
##' \code{ewmeas} and \code{ewcitmeas} are data frames containing weekly
##' reported cases of measles in England and Wales.  \code{ewmeas} records the
##' total measles reports for the whole country, 1948--1966.  One questionable
##' data point has been replaced with an NA.  \code{ewcitmeas} records the
##' incidence in seven English cities 1948--1987.  These data were kindly
##' provided by Ben Bolker, who writes: \dQuote{Most of these data have been
##' manually entered from published records by various people, and are prone to
##' errors at several levels. All data are provided as is; use at your own
##' risk.}
##'
##'
##' @name measles
##' @rdname measles
##' @aliases LondonYorke ewmeas ewcitmeas
##' @docType data
##' @references W. P. London and J. A. Yorke, Recurrent Outbreaks of Measles,
##' Chickenpox and Mumps: I. Seasonal Variation in Contact Rates, American
##' Journal of Epidemiology, 98:453--468, 1973.
##' @keywords datasets
##' @examples
##'
##' plot(cases~time,data=LondonYorke,subset=disease=="measles",type='n',main="measles",bty='l')
##' lines(cases~time,data=LondonYorke,subset=disease=="measles"&town=="Baltimore",col="red")
##' lines(cases~time,data=LondonYorke,subset=disease=="measles"&town=="New York",col="blue")
##' legend("topright",legend=c("Baltimore","New York"),lty=1,col=c("red","blue"),bty='n')
##'
##' plot(
##'      cases~time,
##'      data=LondonYorke,
##'      subset=disease=="chickenpox"&town=="New York",
##'      type='l',col="blue",main="chickenpox, New York",
##'      bty='l'
##'     )
##'
##' plot(
##'      cases~time,
##'      data=LondonYorke,
##'      subset=disease=="mumps"&town=="New York",
##'      type='l',col="blue",main="mumps, New York",
##'      bty='l'
##'     )
##'
##' plot(reports~time,data=ewmeas,type='l')
##'
##' plot(reports~date,data=ewcitmeas,subset=city=="Liverpool",type='l')
##'
NULL


##' Two-dimensional discrete-time Ornstein-Uhlenbeck process
##'
##' \code{ou2} is a \sQuote{pomp} object encoding a bivariate discrete-time
##' Ornstein-Uhlenbeck process.
##'
##' If the state process is \eqn{X(t) = (x_{1}(t),x_{2}(t))}, then \deqn{X(t+1)
##' = \alpha X(t) + \sigma \epsilon(t),} where \eqn{\alpha} and \eqn{\sigma}
##' are 2x2 matrices, \eqn{\sigma} is lower-triangular, and \eqn{\epsilon(t)}
##' is standard bivariate normal.  The observation process is \eqn{Y(t) =
##' (y_1(t),y_2(t))}, where \eqn{y_i(t) \sim \mathrm{normal}(x_i(t),\tau)}.
##' The functions \code{rprocess}, \code{dprocess}, \code{rmeasure},
##' \code{dmeasure}, and \code{skeleton} are implemented using compiled C code
##' for computational speed: see the source code for details.
##'
##' @name ou2
##' @rdname ou2
##' @docType data
##' @seealso \code{\link{pomp}}
##' @keywords models datasets
##' @examples
##'
##' pompExample(ou2)
##' plot(ou2)
##' coef(ou2)
##' x <- simulate(ou2)
##' plot(x)
##' pf <- pfilter(ou2,Np=1000)
##' logLik(pf)
##'
NULL





##' Ricker model with Poisson observations.
##'
##' \code{ricker} is a \sQuote{pomp} object encoding a stochastic Ricker model
##' with Poisson measurement error.
##'
##' The state process is \eqn{N_{t+1} = r N_{t} \exp(-c N_{t}+e_{t})}{N[t+1] =
##' r N[t] exp(-c N[t]+e[t])}, where the \eqn{e_t}{e[t]} are i.i.d. normal
##' random deviates with zero mean and variance \eqn{\sigma^2}{sigma^2}.  The
##' observed variables \eqn{y_t}{y[t]} are distributed as
##' \eqn{\mathrm{Poisson}(\phi N_t)}{Poisson(phi N[t])}.
##'
##' @name ricker
##' @docType data
##' @seealso \code{\link{pomp}}, \code{\link{gompertz}}, and the tutorials on
##' the \href{https://kingaa.github.io/pomp/}{package website}.
##' @keywords datasets models
##' @examples
##'
##' pompExample(ricker)
##' plot(ricker)
##' coef(ricker)
##'
NULL





##' Two-dimensional random-walk process
##'
##' \code{rw2} is a \sQuote{pomp} object encoding a 2-D normal random walk.
##'
##' The random-walk process is fully but noisily observed.
##'
##' @name rw2
##' @docType data
##' @seealso \code{\link{pomp}}, \code{\link{ou2}}
##' @keywords datasets models
##' @examples
##'
##' pompExample(rw2)
##' plot(rw2)
##' x <- simulate(rw2,nsim=10,seed=20348585L,params=c(x1.0=0,x2.0=0,s1=1,s2=3,tau=1))
##' plot(x[[1]])
##'
NULL

##' Compartmental epidemiological models
##'
##' \code{sir} is a \sQuote{pomp} object encoding a simple seasonal SIR model.
##' Simulation is performed using an Euler multinomial approximation.
##'
##' \code{sir2} has the same model implemented using Gillespie's algorithm.
##'
##' \code{bbs} is a nonseasonal SIR model together with data from a 1978
##' outbreak of influenza in a British boarding school.
##'
##' This and similar examples are discussed and constructed in tutorials
##' available on the \href{https://kingaa.github.io/pomp/}{package website}.
##'
##' The boarding school influenza outbreak is described in Anonymous (1978).
##'
##' @name sir_models
##' @rdname sir
##' @aliases sir sir2 bbs
##' @docType data
##' @seealso \code{\link{pomp}} and the tutorials on the
##' \href{https://kingaa.github.io/pomp/}{package website}.
##' @references Anonymous (1978).  Influenza in a boarding school.  British
##' Medical Journal 1:587
##' @keywords datasets models
##' @examples
##'
##' pompExample(sir)
##' plot(sir)
##' plot(simulate(sir))
##' coef(sir)
##'
##' pompExample(sir2)
##' plot(sir2)
##' plot(simulate(window(sir2,end=3)))
##' coef(sir2)
##'
##' pompExample(bbs)
##' plot(bbs)
##' coef(bbs)
##'
NULL
