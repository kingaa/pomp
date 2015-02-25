\name{parus}
\docType{data}
\alias{parus}
\title{POMPs for time series of Parus major population counts from Wytham Wood}
\description{
  \code{pompExample(parus)} constructs a function that constructs \code{pomp} objects implementing each of several population-dynamics models and abundance data from Parus major in Wytham Wood.

  The model has one state variable, \eqn{N_t}, the true abundance.
  There is one observable, \eqn{\mathrm{pop}_t}.

  The function defined is: \preformatted{
  parus(proc = c("Gompertz", "Ricker"),
        meas = c("lognormal","Poisson","negbin"))}
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
pompExample(parus)
p1 <- parus(proc="Ricker",meas="negbin") 
plot(p1)
plot(simulate(p1))

p2 <- parus(proc="Gompertz",meas="lognormal")
plot(simulate(p2))
pf <- pfilter(p2,Np=1000)
}
\seealso{\code{vignette("parus")}}
\keyword{models,datasets}
