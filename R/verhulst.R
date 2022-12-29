##' Verhulst-Pearl model
##'
##' The Verhulst-Pearl (logistic) model of population growth.
##'
##' A stochastic version of the Verhulst-Pearl logistic model.
##' This evolves in continuous time, according to the stochastic differential equation
##' \deqn{dn_t = r\,n_t\,\left(1-\frac{n_t}{K}\right)\,dt+\sigma\,n_t\,dW_t.}{dn[t] = r n[t] (1-n[t]/K) dt + sigma n[t] dW[t].}
##'
##' Numerically, we simulate the stochastic dynamics using an Euler approximation.
##'
##' The measurements are assumed to be log-normally distributed:
##' \deqn{N_t \sim \mathrm{Lognormal}\left(\log{n_t},\tau\right).}{N[t] ~ Lognormal(log(n[t]),tau).}
##'
##' @name verhulst
##' @rdname verhulst
##' @docType data
##' @family pomp examples
##'
##' @param r intrinsic growth rate
##' @param K carrying capacity
##' @param sigma environmental process noise s.d.
##' @param tau measurement error s.d.
##' @param n_0 initial condition
##' @param dt Euler timestep
##'
##' @return
##' A \sQuote{pomp} object containing the model and simulated data.
##' The following basic components are included in the \sQuote{pomp} object:
##' \sQuote{rinit}, \sQuote{rprocess}, \sQuote{rmeasure}, \sQuote{dmeasure}, and \sQuote{skeleton}.
##'
##' @example examples/verhulst.R
##'
NULL

##' @rdname verhulst
##' @export
verhulst <- function (
  n_0 = 10000, K = 10000, r = 0.9,
  sigma = 0.4, tau = 0.1, dt = 0.01)
{

  simulate(
    times=seq(0.1,by=0.1,length=1000),
    t0=0,
    params=c(n_0=n_0,K=K,r=r,sigma=sigma,tau=tau),
    rprocess=euler(
      step.fun=Csnippet("
        n = rnorm(n+r*n*(1-n/K)*dt,sigma*n*sqrt(dt));
      "
      ),
      delta.t=dt
    ),
    emeasure=Csnippet("E_N = n*exp(tau*tau/2);"),
    vmeasure=Csnippet("
      double et = exp(tau*tau);
      V_N_N = n*n*et*(et-1);"),
    rmeasure=Csnippet("N = rlnorm(log(n),tau);"),
    dmeasure=Csnippet("lik = dlnorm(N,log(n),tau,give_log);"),
    skeleton=vectorfield(Csnippet("Dn = r*n*(1-n/K);")),
    rinit=Csnippet("n = n_0;"),
    paramnames=c("r","K","tau","sigma","n_0"),
    statenames=c("n"),
    obsnames="N",
    seed=73658676L
  )

}
