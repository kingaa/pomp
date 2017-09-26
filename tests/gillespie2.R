library(magrittr)
library(pomp)
library(ggplot2)
library(reshape2)

## tests of argument checking

v <- cbind(death = c(-1,1))

works <- function(times = c(1), t0 = 0, mu = 0.001, N_0 = 1) {
    data <- data.frame(time = times, reports = NA)
    rmeasure <- Csnippet("reports = ct;")
    rprocess <- gillespie.sim(Csnippet("rate = mu * N;"), v = v)
    initializer <- function(params, t0, ...) {
        c(N=N_0,ct=12)
    }
    pomp(data = data, times = "time", t0 = t0, params = c(mu=mu),
         rprocess = rprocess, rmeasure = rmeasure, initializer = initializer,
         zeronames="ct", paramnames=c("mu"), statenames=c("N","ct"))
}
mwe <- simulate(works(), seed=1L)

d <- cbind(c(1,0))
try(mwe %>% pomp(rprocess=gillespie.sim(Csnippet("rate = mu * N;"),v=v, d=d),
                 zeronames="ct", paramnames=c("mu"), statenames=c("N","ct")) %>%
    simulate() %>% invisible())

vdup <- cbind(death = c(N=-1,N=1))
try(mwe %>% pomp(rprocess=gillespie.sim(Csnippet("rate = mu * N;"),v=vdup),
                 zeronames="ct", paramnames=c("mu"), statenames=c("N","ct")) %>%
    simulate() %>% invisible())

try(mwe %>% pomp(rprocess=gillespie.sim(Csnippet("rate = mu * N;"),v=vdup),
                 zeronames="ct", paramnames=c("mu"),
                 statenames=c("N","ct", "N")) %>% simulate() %>% invisible())

try(mwe %>% pomp(rprocess=gillespie.sim(Csnippet("rate = mu * N;"),v=vdup),
                 zeronames="ct", paramnames=c("mu", "mu"),
                 statenames=c("N","ct")) %>% simulate() %>% invisible())

try(mwe %>% pomp(rprocess=gillespie.sim(Csnippet("rate = mu * N;"),v=vdup),
                 zeronames=c("ct", "ct"), paramnames=c("mu"),
                 statenames=c("N","ct")) %>% simulate() %>% invisible())

try(mwe %>% pomp(rprocess=gillespie.hl.sim(list(Csnippet("rate = mu * N;"), c(N=-1, ct=1), "bob")),
                 zeronames="ct", paramnames=c("mu"),
                 statenames=c("N","ct")) %>% simulate() %>% invisible())

try(mwe %>% pomp(rprocess=gillespie.hl.sim(.pre=function(x)x,list(Csnippet("rate = mu * N;"), c(N=-1, ct=1))),
             zeronames="ct", paramnames=c("mu"),
             statenames=c("N","ct")) %>% simulate() %>% invisible())

try(mwe %>% pomp(rprocess=gillespie.hl.sim(list(function(N) mu*N, c(N=-1, ct=1))),
                 zeronames="ct", paramnames=c("mu"),
                 statenames=c("N","ct")) %>% simulate() %>% invisible())

mwe %>% pomp(rprocess=gillespie.hl.sim(list(Csnippet("rate = mu * N;"), c(N=-1, ct=1))),
             zeronames="ct", paramnames=c("mu"),
             statenames=c("N","ct")) %>% simulate() %>% invisible()

try(mwe %>% pomp(rprocess=gillespie.hl.sim(list(3L, c(N=-1, ct=1))),
                 zeronames="ct", paramnames=c("mu"),
                 statenames=c("N","ct")) %>% simulate() %>% invisible())

try(mwe %>% pomp(rprocess=gillespie.hl.sim(list(c("rate = mu * N;", "mistake"),
                                                c(N=-1, ct=1))),
                 zeronames="ct", paramnames=c("mu"),
                 statenames=c("N","ct")) %>% simulate() %>% invisible())

try(mwe %>% pomp(rprocess=gillespie.hl.sim(list(c("rate = mu * N;"),
                                                c(N="-1", ct=1))),
                 zeronames="ct", paramnames=c("mu"),
                 statenames=c("N","ct")) %>% simulate() %>% invisible())

try(mwe %>% pomp(rprocess=gillespie.hl.sim(list(c("rate = mu * N;"), c(N=-1))),
                 zeronames="ct", paramnames=c("mu"),
                 statenames=c("N","ct")) %>% simulate() %>% invisible())

try(mwe %>% pomp(rprocess=gillespie.hl.sim(list(c("rate = mu * N;"),
                                                c(N=-1, Ct=1))),
                 zeronames="ct", paramnames=c("mu"),
                 statenames=c("N","ct")) %>% simulate() %>% invisible())

try(mwe %>% pomp(rprocess=gillespie.hl.sim(list(c("rate = mu * N;"),
                                                c(-1, 1))),
                 zeronames="ct", paramnames=c("mu"),
                 statenames=c("N","ct")) %>% simulate() %>% invisible())

## compare simulation output of various model specifications that should all be statistically equivalent

png(filename="gillespie2-%02d.png",res=100)

c(gamma=12,mu=1/70,iota=0.1,
  beta1=330,beta2=410,beta3=490,
  beta.sd=0,
  rho=0.1,
  S_0=0.07,I_0=1e-4,R_0=0.93,
  pop=1000000
) -> params

rate.fun <- function(j, x, t, params, covars, ...) {
  switch(
    j,
    params["mu"]*x["N"],            # birth
    params["mu"]*x["S"],            # susceptible death
    {                               # infection
      beta <- sum(covars*params[c("beta1","beta2","beta3")])
      (beta*x["I"]+params["iota"])*x["S"]/x["N"]
    },
    params["mu"]*x["I"],            # infected death
    params["gamma"]*x["I"],         # recovery
    params["mu"]*x["R"],           # recovered death
    stop("unrecognized event ",j)
  )
}

rate.fun.snip <- Csnippet("
  double beta;
  switch (j) {
  case 1: 			// birth
    rate = mu * N;
    break;
  case 2:			// susceptible death
    rate = mu * S;
    break;
  case 3:			// infection
    beta = seas_1 * beta1 + seas_2 * beta2 + seas_3 * beta3;
    rate = (beta * I + iota) * S / N;
    break;
  case 4:			// infected death
    rate = mu * I;
    break;
  case 5:			// recovery
    rate = gamma * I;
    break;
  case 6:			// recovered death
    rate = mu * R;
    break;
  default:
    error(\"unrecognized event %d\",j);
    break;
  }")

cbind(
  birth=c(1,0,0,1,0),
  sdeath=c(-1,0,0,-1,0),
  infection=c(-1,1,0,0,0),
  ideath=c(0,-1,0,-1,0),
  recovery=c(0,-1,1,0,1),
  rdeath=c(0,0,-1,-1,0)
) -> Vmatrix

data.frame(
  time=seq(from=0,to=2,by=1/52),
  reports=NA
) %>%
  pomp(
    times="time",
    t0=0,
    rprocess=gillespie.sim(rate.fun=rate.fun,v=Vmatrix,hmax=1/52/10),
    zeronames=c("cases"),
    covar=data.frame(
      t=seq(0,2,by=1/52/10),
      seas=periodic.bspline.basis(
        seq(0,2,by=1/52/10),
        degree=3,period=1,nbasis=3)),
    tcovar="t",
    measurement.model=reports~binom(size=cases,prob=rho),
    initializer=function(params, t0, ...){
      comp.names <- c("S","I","R")
      icnames <- paste(comp.names,"0",sep="_")
      snames <- c("S","I","R","N","cases")
      fracs <- params[icnames]
      x0 <- numeric(length(snames))
      names(x0) <- snames
      x0["N"] <- params["pop"]
      x0[comp.names] <- round(params['pop']*fracs/sum(fracs))
      x0
    }
  ) %>%
  simulate(params=params,seed=806867104L) -> gsir

gsir %>%
    pomp(rprocess=gillespie.sim(rate.fun=rate.fun.snip,
                                v=Vmatrix,hmax=1/52/10),
         paramnames = names(params),
         statenames = c("S","I","R", "N", "cases"),
         ) %>%
    simulate(seed=806867104L) -> gsir1

all.equal(as.data.frame(gsir), as.data.frame(gsir1))

hl.args <- list(list("rate = mu * N;",    c(S= 1,I= 0,R= 0,N= 1,cases= 0)),
                list("rate = mu * S;",    c(S=-1,I= 0,R= 0,N=-1,cases= 0)),
                list("rate = mu * I;",    c(S= 0,I=-1,R= 0,N=-1,cases= 0)),
                list("rate = mu * R;",    c(S= 0,I= 0,R=-1,N=-1,cases= 0)),
                list("rate = gamma * I;", c(S= 0,I=-1,R= 1,N= 0,cases= 1)),
                list(paste("double beta;",
                           "beta = seas_1 * beta1 + seas_2 * beta2 + seas_3 * beta3;",
                           "rate = (beta * I + iota) * S / N;", sep = "\n"),
                     c(S=-1,I= 1,R= 0,N= 0,cases= 0)))

gsir %>%
    pomp(rprocess=do.call(gillespie.hl.sim, c(hl.args, list(hmax=1/52/10))),
         paramnames = names(params),
         statenames = c("S","I","R", "N", "cases"),
         ) %>%
    simulate(seed=806867100L) -> gsirhl

hl.args1 <- list(list("rate = mu * N;",    c(R= 0,S= 1,I= 0,cases =0,N= 1)),
                 list("rate = mu * S;",    c(S=-1,I= 0,R= 0,N=-1,cases= 0)),
                 list("rate = mu * R;",    c(S= 0,I= 0,R=-1,N=-1,cases= 0)),
                 list("rate = gamma * I;", c(S= 0,I=-1,R= 1,N= 0,cases= 1)),
                 list("rate = mu * I;",    c(S= 0,I=-1,R= 0,N=-1,cases= 0)),
                list(paste("double beta;",
                           "beta = seas_1 * beta1 + seas_2 * beta2 + seas_3 * beta3;",
                           "rate = (beta * I + iota) * S / N;", sep = "\n"),
                     c(cases= 0,S=-1,I= 1,R= 0,N= 0)))

gsir %>%
    pomp(rprocess=do.call(gillespie.hl.sim, c(hl.args1, list(hmax=1/52/10))),
         paramnames = names(params),
         statenames = c("S","I","R", "N", "cases"),
         ) %>%
    simulate(seed=806867101L) -> gsirhl1

hl.args2 <- list(list("rate = mu * N;",    c(S= 1,I= 0,R= 0,N= 1,cases= 0)),
                list("rate = mu * S;",    c(S=-1,I= 0,R= 0,N=-1,cases= 0)),
                list("rate = mu * I;",    c(S= 0,I=-1,R= 0,N=-1,cases= 0)),
                list("rate = mu * R;",    c(S= 0,I= 0,R=-1,N=-1,cases= 0)),
                list("rate = gamma * I;", c(S= 0,I=-1,R= 1,N= 0,cases= 1)),
                list(paste("beta = seas_1 * beta1 + seas_2 * beta2 + seas_3 * beta3;",
                           "rate = (beta * I + iota) * S / N;", sep = "\n"),
                     c(S=-1,I= 1,R= 0,N= 0,cases= 0)),
                .pre = "double beta;")

gsir %>%
    pomp(rprocess=do.call(gillespie.hl.sim, c(hl.args2, list(hmax=1/52/10))),
         paramnames = names(params),
         statenames = c("S","I","R", "N", "cases"),
         ) %>%
    simulate(seed=806867102L) -> gsirhl2

hl.args3 <- list(list("rate = mu * N / 2;",    c(S= 1,I= 0,R= 0,N= 1,cases= 0)),
                list("rate = mu * S / 2;",    c(S=-1,I= 0,R= 0,N=-1,cases= 0)),
                list("rate = mu * I/ 2;",    c(S= 0,I=-1,R= 0,N=-1,cases= 0)),
                list("rate = mu * R/ 2;",    c(S= 0,I= 0,R=-1,N=-1,cases= 0)),
                list("rate = gamma * I / 2;", c(S= 0,I=-1,R= 1,N= 0,cases= 1)),
                list(paste("beta = seas_1 * beta1 + seas_2 * beta2 + seas_3 * beta3;",
                           "rate = (beta * I + iota) * S / N / 2;", sep = "\n"),
                     c(S=-1,I= 1,R= 0,N= 0,cases= 0)),
                .pre = "double beta;", .post = "rate *= 2;")

gsir %>%
    pomp(rprocess=do.call(gillespie.hl.sim, c(hl.args3, list(hmax=1/52/10))),
         paramnames = names(params),
         statenames = c("N", "S","I","R", "cases"),
         ) %>%
    simulate(seed=806867103L) -> gsirhl3

list(gill.sim.R=as.data.frame(gsir),
     gill.hl.sim=as.data.frame(gsirhl),
     gill.hl.sim1=as.data.frame(gsirhl1),
     gill.hl.sim2=as.data.frame(gsirhl2),
     gill.hl.sim3=as.data.frame(gsirhl3)) %>%
  melt(id="time") %>%
  subset(variable=="reports") %>%
  ggplot(aes(x=time,y=value,color=L1))+
  labs(color="",y="reports",title="various implementations of same SIR model")+
  geom_line()+
  theme_bw()+theme(legend.position=c(0.8,0.8))
