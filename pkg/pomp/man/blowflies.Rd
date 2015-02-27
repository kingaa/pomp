\name{blowflies}
\alias{blowflies}
\alias{blowflies1}
\alias{blowflies2}
\docType{data}
\title{Model for Nicholson's blowflies.}
\description{
  \code{blowflies1} and \code{blowflies2} are \code{pomp} objects encoding stochastic delay-difference models.
}
\examples{
pompExample(blowflies)
plot(blowflies1)
plot(blowflies2)
}
\details{
  The data are from "population I", a control culture in one of A. J. Nicholson's experiments with the Australian sheep-blowfly \emph{Lucilia cuprina}.
  The experiment is described on pp. 163--4 of Nicholson (1957).
  Unlimited quantities of larval food were provided;
  the adult food supply (ground liver) was constant at 0.4g per day.
  The data were taken from the table provided by Brillinger et al. (1980).

  The models are discrete delay equations:
  \deqn{R(t+1) \sim \mathrm{Poisson}(P N(t-\tau) \exp{(-N(t-\tau)/N_{0})} e(t+1) {\Delta}t)}{%
    R[t+1] ~ Poisson(P N[t-tau] exp(-N[t-tau]/N0) e[t+1] dt)}
  \deqn{S(t+1) \sim \mathrm{binomial}(N(t),\exp{(-\delta \epsilon(t+1) {\Delta}t)})}{%
    S[t+1] ~ binomial(N[t],exp(-delta eps[t+1] dt))}
  \deqn{N(t) = R(t)+S(t)}{N[t]=R[t]+S[t]}
  where \eqn{e(t)}{e[t]} and \eqn{\epsilon(t)}{eps[t]} are Gamma-distributed i.i.d. random variables with mean 1 and variances \eqn{{\sigma_p^2}/{{\Delta}t}}{sigma.p^2/dt}, \eqn{{\sigma_d^2}/{{\Delta}t}}{sigma.d^2/dt}, respectively.
  \code{blowflies1} has a timestep (\eqn{{\Delta}t}{dt}) of 1 day, and \code{blowflies2} has a timestep of 2 days.
  The process model in \code{blowflies1} thus corresponds exactly to that studied by Wood (2010).
  The measurement model in both cases is taken to be
  \deqn{y(t) \sim \mathrm{negbin}(N(t),1/\sigma_y^2)}{y[t] ~ negbin(N[t],1/sigma.y^2)}, i.e., the observations are assumed to be negative-binomially distributed with mean \eqn{N(t)}{N[t]} and variance \eqn{N(t)+(\sigma_y N(t))^2}{N[t]+(sigma.y N[t])^2}.
  
  Do

  \code{file.show(system.file("examples","blowflies.R",package="pomp"))}

  to view the code that constructs these pomp objects.
}
\seealso{\code{\link{pomp}}}
\references{
  A. J. Nicholson (1957)
  The self-adjustment of populations to change.
  Cold Spring Harbor Symposia on Quantitative Biology, \bold{22}, 153--173.
  
  Y. Xia and H. Tong (2011)
  Feature Matching in Time Series Modeling.
  \emph{Statistical Science} \bold{26}, 21--46.

  E. L. Ionides (2011)
  Discussion of ``Feature Matching in Time Series Modeling'' by Y. Xia and H. Tong.
  \emph{Statistical Science} \bold{26}, 49--52.

  S. N. Wood (2010)
  Statistical inference for noisy nonlinear ecological dynamic systems.
  \emph{Nature} \bold{466}, 1102--1104.
 
  W. S. C. Gurney, S. P. Blythe, and R. M. Nisbet (1980)
  Nicholson's blowflies revisited.
  \emph{Nature} \bold{287}, 17--21.

  D. R. Brillinger, J. Guckenheimer, P. Guttorp and G. Oster (1980)
  Empirical modelling of population time series: The case of age and density dependent rates.
  in G. Oster (ed.), Some Questions in Mathematical Biology, vol. 13, pp. 65--90.
  American Mathematical Society, Providence.
}
\keyword{models}
\keyword{datasets}
