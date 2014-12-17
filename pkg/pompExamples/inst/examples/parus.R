require(pomp)

dat <- 'year,pop
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
1986,211
'

dat <- read.csv(text=dat)

parus.example <- function (proc = c("Gompertz","Ricker"),
                           meas = c("lognormal","Poisson","negbin")) {

  proc <- match.arg(proc)
  meas <- match.arg(meas)

  pomp(
       data=dat,
       times="year",
       t0=1960,
       params=c(K=190,r=2.7,sigma=0.2,theta=0.05,N.0=148),
       rprocess=discrete.time.sim(
         step.fun=switch(proc,
           Gompertz="_parus_gompertz_simulator",
           Ricker="_parus_ricker_simulator",
           stop("unrecognized value of ",sQuote("proc"))
           ),
         delta.t=1
         ),
       skeleton=switch(proc,
         Gompertz="_parus_gompertz_skeleton",
         Ricker="_parus_ricker_skeleton",
         stop("unrecognized value of ",sQuote("proc"))
         ),
       skeleton.type="map",
       skelmap.delta.t=1,
       rmeasure=switch(meas,
         lognormal="_parus_lognormal_rmeasure",
         Poisson="_parus_poisson_rmeasure",
         negbin="_parus_nbinom_rmeasure",
         stop("unrecognized value of ",sQuote("meas"))
         ),
       dmeasure=switch(meas,
         lognormal="_parus_lognormal_dmeasure",
         Poisson="_parus_poisson_dmeasure",
         negbin="_parus_nbinom_dmeasure",
         stop("unrecognized value of ",sQuote("meas"))
         ),
       paramnames=c("r","K","sigma","theta"),
       statenames=c("N"),
       obsnames=c("pop"),
       parameter.transform=function(params,...){
         exp(params)
       },
       parameter.inv.transform=function(params,...){
         log(params)
       },
       PACKAGE="pompExamples"
       )
}

parus <- parus.example(proc=proc,meas=meas)

c("parus")
