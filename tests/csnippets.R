library(pomp)
set.seed(588400992L)

library(ggplot2)
library(plyr)
library(reshape2)
library(magrittr)
theme_set(theme_bw())

parus.dat <- read.csv(text="
                      year,P
                      1960,148
                      1961,258
                      1962,185
                      1963,170
                      1964,267
                      1965,239
                      1966,196
                      1967,132
                      1968,167
                      1969,186
                      1970,128
                      1971,227
                      1972,174
                      1973,177
                      1974,137
                      1975,172
                      1976,119
                      1977,226
                      1978,166
                      1979,161
                      1980,199
                      1981,306
                      1982,206
                      1983,350
                      1984,214
                      1985,175
                      1986,211"
                      )

step.fun <- Csnippet("
  double dW = rnorm(0,sqrt(dt));
  N += r*N*(1-N/K)*dt+sigma*N*dW;
")

rmeas <- Csnippet("
  P = rpois(N);
")

dmeas <- Csnippet("
  lik = dpois(P,N,give_log);
")

parus <- pomp(data=parus.dat,time="year",t0=1959,
              rprocess=euler.sim(step.fun=step.fun,delta.t=1/365),
	      dmeasure=dmeas,rmeasure=rmeas,
              statenames="N",paramnames=c("r","K","sigma"),
	      cdir=getwd(),cfile="parus")

sim <- simulate(parus,params=c(r=0.2,K=200,sigma=0.5,N.0=200),
                nsim=10,obs=TRUE,states=TRUE)

pf <- pfilter(parus,Np=1000,params=c(r=0.2,K=200,sigma=0.5,N.0=200))
logLik(pf)
