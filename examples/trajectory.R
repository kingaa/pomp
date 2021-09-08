\donttest{
  ## The basic components needed to compute trajectories
  ## of a deterministic dynamical system are
  ## rinit and skeleton.

  ## The following specifies these for a simple continuous-time
  ## model: dx/dt = r (1+e cos(t)) x

  trajectory(
    t0 = 0, times = seq(1,30,by=0.1),
    rinit = function (x0, ...) {
      c(x = x0)
    },
    skeleton = vectorfield(
      function (r, e, t, x, ...) {
        c(x=r*(1+e*cos(t))*x)
      }
    ),
    params = c(r=1,e=3,x0=1)
  ) -> po

  plot(po,log='y')

  ## In the case of a discrete-time skeleton,
  ## we use the 'map' function.  For example,
  ## the following computes a trajectory from
  ## the dynamical system with skeleton
  ## x -> x exp(r sin(omega t)).

  trajectory(
    t0 = 0, times=seq(1,100),
    rinit = function (x0, ...) {
      c(x = x0)
    },
    skeleton = map(
      function (r, t, x, omega, ...) {
        c(x=x*exp(r*sin(omega*t)))
      },
      delta.t=1
    ),
    params = c(r=1,x0=1,omega=4)
  ) -> po

  plot(po)

}
