## blowfly model, with general dt
## here, set up for dt=1 and dt=2
## dt is hard-coded, and initial values are customized for each dt

library(pomp2)

## following xia and tong, the delay is treated as fixed at 14 days
## xia and tong claim to be using tau=8 bidays, but on inspection
## their Euler method is really tau=7 bidays


raw.data <- subset(blowflies,set==4,select=-set)
raw.data$y <- raw.data$count
raw.data$count <- NULL

pomp(
  data=subset(raw.data,day>14&day<400),
  times="day",
  t0=14,
  cfile="blowfly1_source",
  compile=FALSE,
  rprocess=discrete_time(
    step.fun=Csnippet("
      double *N = &N1;
      int tau = 14, k;
      e = rgammawn(sigma_P,dt)/dt;
      eps = rgammawn(sigma_d,dt)/dt;
      R = rpois(P*N[tau]*exp(-N[tau]/N0)*dt*e);
      S = rbinom(N[0],exp(-delta*dt*eps));
      for (k = tau; k > 0; k--) N[k] = N[k-1];
      N[0] = R+S;"
    ),
    delta.t=1
  ),
  dmeasure=Csnippet("
      double size = 1.0/sigma_y/sigma_y;
      lik = dnbinom_mu(y,size,N1,give_log);"
  ),
  rmeasure=Csnippet("
      double size = 1.0/sigma_y/sigma_y;
      y = rnbinom_mu(size,N1);"
  ),
  partrans=parameter_trans(
    log=c("P","delta","N0","sigma.P","sigma.d","sigma.y")
  ),
  paramnames=c("P","N0","delta","sigma.P","sigma.d","sigma.y"),
  statenames=c("N1","R","S","e","eps"),
  y.init=with( ## initial data
    raw.data,
    approx(x=day,y=y,xout=seq(from=0,to=14,by=1),rule=2)$y
  ),
  #     y.init=c(948, 948, 942, 930, 911, 885, 858, 833.7, 801, 748.3, 676, 589.8, 504, 434.9, 397),
  rinit=function (params, t0, y.init, ...) {
    ntau <- length(y.init)
    n <- y.init[ntau:1]
    names(n) <- paste("N",seq_len(ntau),sep="")
    c(n,R=0,S=0,e=0,eps=0)
  }
) -> blowflies1

pomp(
  blowflies1,
  cfile="blowfly2_source",
  compile=FALSE,
  rprocess=discrete_time(
    step.fun=Csnippet("
      double *N = &N1;
      int tau = 7, k;
      e = rgammawn(sigma_P,dt)/dt;
      eps = rgammawn(sigma_d,dt)/dt;
      R = rpois(P*N[tau]*exp(-N[tau]/N0)*dt*e);
      S = rbinom(N[0],exp(-delta*dt*eps));
      for (k = tau; k > 0; k--) N[k] = N[k-1];
      N[0] = R+S;"
    ),
    delta.t=2
  ),
  y.init=with( ## initial data
    raw.data,
    approx(x=day,y=y,xout=seq(from=0,to=14,by=2),rule=2)$y
  ),
  #y.init=c(948, 942, 911, 858, 801, 676, 504, 397),
  paramnames=c("P","N0","delta","sigma.P","sigma.d","sigma.y"),
  statenames=c("N1","R","S","e","eps"),
) -> blowflies2

## mle from search to date
coef(blowflies1,transform=TRUE) <- c(
  P = 1.189,
  delta = -1.828,
  N0 = 6.522,
  sigma.P = 0.301,
  sigma.d = -0.292,
  sigma.y = -3.625
)

## mle from search to date
coef(blowflies2,transform=TRUE) <- c(
  P = 1.005,
  delta = -1.75,
  N0 = 6.685,
  sigma.P = 0.366,
  sigma.d = -0.274,
  sigma.y = -4.524
)

c("blowflies1","blowflies2")
