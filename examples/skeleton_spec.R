\donttest{
  ## Starting with an existing pomp object,
  ## e.g., the continuous-time Verhulst-Pearl model,

  verhulst() -> po
  
  ## we add or change the deterministic skeleton
  ## using the 'skeleton' argument in any 'pomp'
  ## estimation function, or in the 'pomp' constructor
  ## itself). Here, we pass the skeleton specification
  ## to 'pomp' as an R function.
  ## Since this is a continuous-time POMP, the
  ## skeleton is a vectorfield.

  po %>%
    pomp(
      skeleton=vectorfield(
        function(r, K, n, ...) {
          c(n=r*n*(1-n/K))
        }
      )
    ) %>%
    trajectory(format="data.frame") -> traj

  ## We can also pass it as a C snippet:

  po %>%
    traj_objfun(
      skeleton=vectorfield(Csnippet("Dn=r*n*(1-n/K);")),
      paramnames=c("r","K"),
      statenames="n"
    ) -> ofun

  ofun()

  ## For a discrete-time POMP, the deterministic skeleton
  ## is a map.  For example,

  gompertz() -> po

  po %>%
    traj_objfun(
      skeleton=map(
        Csnippet("
          double dt = 1.0;
          double s = exp(-r*dt);
          DX = pow(K,(1-s))*pow(X,s);"
        ), delta.t=1
      ),
      paramnames=c("r","K"),
      statenames=c("X")
    ) -> ofun

  ofun()

}
