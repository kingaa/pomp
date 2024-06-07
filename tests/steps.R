options(digits=3)

suppressPackageStartupMessages({
  library(tidyr)
  library(dplyr)
})
library(pomp)

tm <- freeze(sort(runif(n=20,max=50)),seed=1930502785L)
simulate(
  times=tm,
  t0=0,
  rprocess=discrete_time(
    Csnippet("
      double dW = rnorm(0,sqrt(dt));
      N += r*N*(1-N)*dt+sigma*dW;
      e += dW;
      step += 1;
    "),
    delta.t=0.1
  ),
  skeleton=map(
    Csnippet("
      double dt = 0.1;
      DN = N+r*N*(1-N)*dt;
      De = 0;
      Dstep = step+1;
    "),
    delta.t=0.1
  ),
  rinit=Csnippet("e=0; N = N_0; step = 0;"),
  accumvars="step",
  params=c(r=0.5,N_0=0.5,sigma=0.1,phi=10,c=1),
  paramnames=c("sigma","r","N_0"),
  statenames=c("N","e","step")
) -> sm
sm |> trajectory() -> tj
stopifnot(
  all.equal(
    diff(floor(time(sm,t0=TRUE)/0.1)),
    as.numeric(states(sm,"step"))
  ),
  all.equal(
    floor(diff(time(tj,t0=TRUE))/0.1),
    as.numeric(states(tj,"step"))
  )
)

tm <- freeze(sort(runif(n=100,max=2)),seed=342643819)
sir() |>
  simulate(
    times=tm,
    rprocess=euler(Csnippet("
      double rate[6];		  // transition rates
      double trans[6];		// transition numbers
      double beta;
      double dW;
      beta = dot_product(3,&beta1,&seas_1);
      dW = rgammawn(beta_sd,dt);
      rate[0] = mu*pop;
      rate[1] = (iota+beta*I*dW/dt)/pop;
      rate[2] = mu;
      rate[3] = gamma;
      rate[4] = mu;
      rate[5] = mu;
      trans[0] = rpois(rate[0]*dt);	// births are Poisson
      reulermultinom(2,S,&rate[1],dt,&trans[1]);
      reulermultinom(2,I,&rate[3],dt,&trans[3]);
      reulermultinom(1,R,&rate[5],dt,&trans[5]);
      S += trans[0]-trans[1]-trans[2];
      I += trans[1]-trans[3]-trans[4];
      R += trans[3]-trans[5];
      cases += trans[3];
      steps += 1;
      "),
      delta.t=0.01
    ),
    rinit=Csnippet("
      double m = pop/(S_0+I_0+R_0);
      S = nearbyint(m*S_0);
      I = nearbyint(m*I_0);
      R = nearbyint(m*R_0);
      cases = 0;
      W = 0;
      steps = 0;
      "),
    statenames=c("S","I","R","cases","steps","W"),
    accumvars=c("cases","steps"),
    paramnames=c(
      "mu","gamma","iota","beta1","beta_sd","pop",
      "S_0","I_0","R_0"
    ),
    seed=1232934371
  ) -> sm

stopifnot(
  all.equal(
    ceiling(diff(time(sm,t0=TRUE))/0.01),
    as.numeric(states(sm,"steps"))
  )
)

simulate(
  t0=0,
  times=c(0,0,0,1,3,3,7),
  rinit=\(t0,...)c(x=t0,n=0),
  rprocess=onestep(\(t,x,n,delta.t,...) c(x=t+delta.t,n=n+1)),
  rmeasure=\(x,n,...) c(y1=x,y2=n),
  dmeasure=\(...,log) if (log) 0 else 1
) -> s1

simulate(
  t0=0,
  times=c(0,0,0,1,3,3,7),
  rinit=Csnippet("x=t; n=0;"),
  rprocess=onestep(Csnippet("x=t+dt; n=n+1;")),
  rmeasure=Csnippet("y1 = x; y2 = n;"),
  dmeasure=Csnippet("lik = (give_log) ? 0 : 1;"),
  statenames=c("x","n"),
  obsnames=c("y1","y2")
) -> s2

s1 |>
  as.data.frame() |>
  with(
    stopifnot(
      time==x,
      n==seq_along(time)
    )
  )

s2 |>
  pfilter(Np=1,save.states=TRUE) |>
  saved_states(format="list") |>
  melt() |>
  pivot_wider() |>
  mutate(time=time(s2)[.L1]) |>
  with(
    stopifnot(
      time==x,
      n==seq_along(time)
    )
  )
