library(pomp)
set.seed(807969746L)

gompertz() -> gompertz
theta <- coef(gompertz)
theta["X_0"] <- 2
gompertz |>
  simulate(
    dinit=function(X, X_0, t0, ...) {
      if (X==X_0) 0 else -Inf
    }
  ) -> sm

sm |> dinit(x=rinit(sm,nsim=5))
sm |> dinit(params=theta,x=rinit(sm),t0=3)
sm |> dinit(params=cbind(parmat(theta,3),coef(sm)),x=rinit(sm))
sm |> dinit(params=theta,x=rinit(sm,nsim=5))
try(sm |> dinit(params=cbind(parmat(theta,2),coef(sm)),x=rinit(sm,nsim=2)))
try(sm |> dinit(params=cbind(parmat(theta,3),coef(sm)),x=rinit(sm,nsim=5)))
sm |> dinit(params=cbind(parmat(theta,3),coef(sm)),x=rinit(sm,nsim=8))

gompertz |>
  simulate(
    rinit=Csnippet("X = rexp(1/X_0);"),
    dinit=Csnippet("loglik = dexp(X,1/X_0,1);"),
    statenames="X",
    paramnames="X_0"
  ) -> sm
sm |> dinit(x=rinit(sm,nsim=5))
sm |> dinit(params=theta,x=rinit(sm),t0=3)
sm |> dinit(params=cbind(parmat(theta,3),coef(sm)),x=rinit(sm))
sm |> dinit(params=theta,x=rinit(sm,nsim=5))
try(sm |> dinit(params=cbind(parmat(theta,2),coef(sm)),x=rinit(sm,nsim=2)))
try(sm |> dinit(params=cbind(parmat(theta,3),coef(sm)),x=rinit(sm,nsim=5)))
sm |> dinit(params=cbind(parmat(theta,3),coef(sm)),x=rinit(sm,nsim=8))

try(dinit(x=rinit(sm)))
try(dinit("bob",x=rinit(sm)))

sir() -> sir
theta <- coef(sir)
sir |>
  simulate(
    dinit=function (S, I, R, S_0, I_0, R_0, t0,
      seas_1, seas_2, seas_3, pop, ...) {
      print(c(t0,seas_1+seas_2+seas_3))
      frac <- c(S,I,R)/(S+I+R)
      if (all(frac==c(S_0,I_0,R_0))) {
        0
      } else {
        -Inf
      }
    }
  ) |>
  dinit(x=rinit(sir),params=parmat(theta,2),log=FALSE)

sir |>
  dinit(x=rinit(sir),params=parmat(theta,2),log=FALSE)

ricker() |>
  dinit(
    x=cbind(c(N=7,e=0),c(N=7,e=1),c(N=4,e=0))
  )
