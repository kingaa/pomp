options(digits=3)
png(filename="steps-%02d.png",res=100)

library(pomp2)
library(magrittr)

ricker() -> ricker

set.seed(54588699L)

po <- ricker
coef(po,"sigma") <- 0
tm <- sort(runif(n=20,max=3))
x <- trajectory(po,times=tm)["N",,]
y <- simulate(po,times=tm,format="arrays")$states["N",,]
stopifnot(identical(x,y))

po <- simulate(ricker,
  times=sort(runif(n=20,max=50)),
  rprocess=discrete_time(
    Csnippet("
      double dW = rnorm(0,sqrt(dt));
      N += r*N*(1-N)*dt+sigma*dW;
      e += dW;
      step += 1;
    "),
    delta.t=1),
  skeleton=NULL,
  rinit=Csnippet("e=0; N = N_0; step = 0;"),
  zeronames="step",
  params=c(r=0.5,N_0=0.5,sigma=0.1,phi=10,c=1),
  paramnames=c("sigma","r","N_0"),statenames=c("N","e","step"))
plot(po)
states(po,"step") %>% table()

sir() -> sir
tm <- sort(runif(n=100,max=2))
sir %>%
  simulate(times=tm,
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
      "),delta.t=0.01),
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
    zeronames=c("cases","steps"),
    paramnames=c("mu","gamma","iota","beta1","beta_sd","pop",
      "S_0","I_0","R_0"),
    seed=1232934371,
    format="arrays"
  ) %>%
  extract2("states") %>%
  extract("steps",,) -> x
plot(diff(c(0,tm)),x)
abline(0,1/0.01)
x %>% table()

dev.off()
