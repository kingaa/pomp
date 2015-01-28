\name{parus}
\docType{data}
\alias{parus}
\title{POMPs for time series of Parus major population counts from Wytham Wood}
\description{
  \code{pompExample(parus,proc,meas)} constructs a \code{pomp} object containing one of several population-dynamics models and abundance data from Parus major in Wytham Wood.

  The model has one state variables, \eqn{N_t}, the true abundance.
  There is one observable, \eqn{\mathrm{pop}_t}.
}
\section{State process}{
  If \code{proc="Gompertz"}, the state process is the stochastic Gompertz model
  \deqn{N_{t} = K^{1-S}\,N_{t-1}^S\,\epsilon_t,}
  where \eqn{S=e^r} and \eqn{\epsilon_t\sim\mathrm{Lognormal}(0,\sigma)} i.i.d.

  If \code{proc="Ricker"}, the state process is the stochastic Ricker model
  \deqn{N_{t} = N_{t-1}\,\exp\left(r\,\left(1-\frac{N_{t-1}}{K}\right)\right)\,\epsilon_t,}
  where, again, \eqn{\epsilon_t\sim\mathrm{Lognormal}(0,\sigma)} i.i.d.
}  
\section{Measurement process}{
  There are three alternative measurement models.

  If \code{meas="Poisson"}, the measurement process is
  \deqn{\mathrm{pop}_{t} \sim \mathrm{Poisson}(N_t).}

  If \code{meas="negbin"}, the measurement process is
  \deqn{\mathrm{pop}_{t} \sim \mathrm{Negbin}(N_t,\theta),}
  i.e., \eqn{\mathrm{pop}_t} has mean \eqn{N_t} and variance \eqn{N_t+\frac{N_t^2}{\theta}}.

  If \code{meas="lognormal"}, the measurement process is
  \deqn{\mathrm{pop}_{t} \sim \mathrm{Lognormal}(\log(N_t),\theta).}
}
\examples{
pompExample(parus,proc="Ricker",meas="negbin")
plot(parus)
plot(simulate(parus))

pompExample(parus,proc="Gompertz",meas="lognormal")
plot(simulate(parus))
}
\seealso{\code{vignette("parus")}}
\keyword{models,datasets}
