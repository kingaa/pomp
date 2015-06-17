## data read from graph in Anonymous (1978) 'Influenza in a boarding school' Brit. Med. J. 1:578.
## cases are recorded with error of +/- 1 case
## 763 boys were at risk, 512 boys spent time away from class
flu <- read.csv(text="
date,confined,convalescent
1978-01-22,1,0
1978-01-23,6,0
1978-01-24,26,0
1978-01-25,73,1
1978-01-26,222,8
1978-01-27,293,16
1978-01-28,258,99
1978-01-29,236,160
1978-01-30,191,173
1978-01-31,124,162
1978-02-01,69,150
1978-02-02,26,89
1978-02-03,11,44
1978-02-04,4,22
",colClasses=c(date='Date'))
flu$day <- flu$date-min(flu$date)+1
units(flu$day) <- "days"
flu$day <- as.numeric(flu$day)

partrans <- "
 TBeta = exp(Beta);
 Tinf_pd = exp(inf_pd);
 Trho = expit(rho);
 Tsfrac = expit(sfrac);
"

paruntrans <- "
 TBeta = log(Beta);
 Tinf_pd = log(inf_pd);
 Trho = logit(rho);
 Tsfrac = logit(sfrac);
"

dmeas <- "
  lik = dpois(confined,rho*R+1e-6,give_log);
"

rmeas <- "
  confined = rpois(rho*R+1e-6);
  convalescent = rpois(rho*C);
"

stochsim <- "
  double t1 = rbinom(S,1-exp(-Beta*I*dt));
  double t2 = rbinom(I,1-exp(-dt/inf_pd));
  double t3 = rbinom(R,1-exp(-dt/conf_pd));
  double t4 = rbinom(C,1-exp(-dt/conv_pd));
  S -= t1;
  I += t1 - t2;
  R += t2 - t3;
  C += t3 - t4;
"

skel <- "
  double dt = 1.0/24.0;
  double t1 = S*(1-exp(-Beta*I*dt));
  double t2 = I*(1-exp(-dt/inf_pd));
  double t3 = R*(1-exp(-dt/conf_pd));
  double t4 = C*(1-exp(-dt/conv_pd));
  DS = S - t1;
  DI = I + t1 - t2;
  DR = R + t2 - t3;
  DC = C + t3 - t4;
"

pomp(
     data=flu[c("day","confined","convalescent")],
     times="day",
     t0=0,
     params=c(
       Beta=0.004,
       inf.pd=0.7,
       conf.pd=sum(flu$confined)/512,
       conv.pd=sum(flu$convalescent)/512,
       rho=0.9,
       sfrac=762/763
       ),
     rprocess=euler.sim(
       step.fun=Csnippet(stochsim),
       delta.t=1/24
       ),
     skeleton=Csnippet(skel),
     skelmap.delta.t=1/24,
     skeleton.type="map",
     rmeasure=Csnippet(rmeas),
     dmeasure=Csnippet(dmeas),
     fromEstimationScale=Csnippet(partrans),
     toEstimationScale=Csnippet(paruntrans),
     obsnames = c("confined","convalescent"),
     statenames=c("S","I","R","C"),
     paramnames=c(
       "Beta",
       "inf.pd","conf.pd","conv.pd",
       "rho","sfrac"
       ),
     initializer=function(params, t0, ...) {
       x0 <- setNames(numeric(4),c("S","I","R","C"))
       S.0 <- round(763*params["sfrac"])
       x0[c("S","I")] <- c(S.0,763-S.0)
       x0
     }
     ) -> bsflu

c("bsflu")
