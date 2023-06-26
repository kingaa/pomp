\donttest{
  ## The familiar Ricker example.
  ## Suppose that for some reason we wish to pass 'phi'
  ## via the userdata facility instead of as a parameter.

  ## C snippet approach:

  simulate(times=1:100,t0=0,
    phi=as.double(100),
    params=c(r=3.8,sigma=0.3,N.0=7),
    rprocess=discrete_time(
      step.fun=Csnippet(r"{
      double e = (sigma > 0.0) ? rnorm(0,sigma) : 0.0;
      N = r*N*exp(-N+e);}"
      ),
      delta.t=1
    ),
    rmeasure=Csnippet(r"{
       double phi = *get_userdata_double("phi");
       y = rpois(phi*N);}"
    ),
    paramnames=c("r","sigma"),
    statenames="N",
    obsnames="y"
  ) -> rick1

  ## The same problem solved using 'globals':
  simulate(times=1:100,t0=0,
    globals=Csnippet("static double phi = 100;"),
    params=c(r=3.8,sigma=0.3,N.0=7),
    rprocess=discrete_time(
      step.fun=Csnippet(r"{
      double e = (sigma > 0.0) ? rnorm(0,sigma) : 0.0;
      N = r*N*exp(-N+e);}"
      ),
      delta.t=1
    ),
    rmeasure=Csnippet("
       y = rpois(phi*N);"
    ),
    paramnames=c("r","sigma"),
    statenames="N",
    obsnames="y"
  ) -> rick2

  ## Finally, the R function approach:

  simulate(times=1:100,t0=0,
    phi=100,
    params=c(r=3.8,sigma=0.3,N_0=7),
    rprocess=discrete_time(
      step.fun=function (r, N, sigma, ...) {
        e <- rnorm(n=1,mean=0,sd=sigma)
        c(N=r*N*exp(-N+e))
      },
      delta.t=1
    ),
    rmeasure=function (phi, N, ...) {
      c(y=rpois(n=1,lambda=phi*N))
    }
  ) -> rick3

}
