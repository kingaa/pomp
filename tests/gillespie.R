library(pomp)
library(magrittr)
library(reshape2)
library(ggplot2)

set.seed(754646834L)

png(filename="gillespie-%02d.png",res=100)

c(gamma=24,mu=1/70,iota=0.1,
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
  const double *BETA = &beta1;
  const double *s = &seas_1;
  int nbasis = *get_pomp_userdata_int(\"nbasis\");
  int k;

  switch (j) {
  case 1:                       // birth
    rate = mu*pop;
    break;
  case 2:                       // susceptible death
    rate = mu*S;
    break;
  case 3:                       // infection
    for (k = 0, beta = 0; k < nbasis; k++) beta += s[k]*BETA[k];
    rate = (beta*I+iota)*S/pop;
    break;
  case 4:                       // infected death
    rate = mu*I;
    break;
  case 5:                       // recovery
    rate = gamma*I;
    break;
  case 6:                       // recovered death
    rate = mu*R;
    break;
  }"
)

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
    initializer=function (params, t0, ...) {
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
  pomp(zeronames=NULL) %>%
  simulate() %>%
  plot(main="Gillespie SIR, no zeroing")

gsir %>%
  pomp(
    rprocess=gillespie.sim(
      rate.fun=rate.fun.snip,
      v=Vmatrix,hmax=1/52/10
    ),
    nbasis=3L,
    paramnames=c("gamma","mu","iota","beta1","beta2","beta3","beta.sd",
      "pop","rho"),
    statenames=c("S","I","R","N","cases")
  ) %>%
  simulate() %>%
  plot(main="Gillespie SIR, with zeroing")

gsir %>%
  pomp(rprocess=gillespie.sim(rate.fun=rate.fun.snip,
    v=Vmatrix,hmax=1/52/10),
    nbasis=3L,
    paramnames=names(params),
    statenames=c("S","I","R","N","cases"),
  ) %>%
  simulate(seed=806867104L) -> gsir1

pompExample(gillespie.sir)
gsir2 <- simulate(gillespie.sir,params=coef(gsir),
  times=time(gsir),t0=timezero(gsir),seed=806867104L)

list(R=as.data.frame(gsir),
  Csnippet=as.data.frame(gsir1),
  C=as.data.frame(gsir2)) %>%
  melt(id="time") %>%
  subset(variable=="reports") %>%
  ggplot(aes(x=time,y=value,color=L1))+
  labs(color="",y="reports",title="C vs R implementations")+
  geom_line()+
  theme_bw()+theme(legend.position=c(0.2,0.8))

try(gsir %>% pomp(rprocess=gillespie.sim(rate.fun=rate.fun,v=as.numeric(Vmatrix))))

stopifnot(
  gsir %>%
    simulate(params=c(gamma=0,mu=0,iota=0,beta1=0,beta2=0,beta3=0,beta.sd=0,
      rho=0.1,S_0=0.07,I_0=1e-4,R_0=0.93,pop=1e6)) %>%
    states() %>% apply(1,diff) %>% equals(0) %>% all())

rate.fun.bad <- function(j, x, t, params, covars, ...) {
  if (t>1) {
    as.numeric(0)
  } else {
    rate.fun(j,x,t,params,covars,...)
  }
}

pomp(gsir,rprocess=gillespie.sim(rate.fun=rate.fun.bad,v=Vmatrix)) %>%
  simulate() %>%
  plot(main="freeze at time 1")

rate.fun.bad <- function(j, x, t, params, covars, ...) {
  if (t>1) {
    -rate.fun(j,x,t,params,covars,...)
  } else {
    rate.fun(j,x,t,params,covars,...)
  }
}

try(pomp(gsir,rprocess=gillespie.sim(rate.fun=rate.fun.bad,v=Vmatrix)) %>% simulate())

rate.fun.bad <- function(j, x, t, params, covars, ...) -1

try(pomp(gsir,rprocess=gillespie.sim(rate.fun=rate.fun.bad,v=Vmatrix)) %>% simulate())

rate.fun.bad <- function(j, x, t, params, covars, ...) c(1,1)

try(pomp(gsir,rprocess=gillespie.sim(rate.fun=rate.fun.bad,v=Vmatrix)) %>% simulate())

try(pomp(gsir,rprocess=gillespie.sim(rate.fun=3,v=Vmatrix)))

create_example <- function(times = c(1,2), t0 = 0, mu = 0.001, N_0 = 1) {
  data <- data.frame(time = times, reports = NA)
  mmodel <- reports ~ pois(ct)
  rate.fun <- function(j, x, t, params, covars, ...) {
    switch(j, params["mu"]*x["N"], stop("unrecognized event ",j))
  }
  rprocess <- gillespie.sim(rate.fun = rate.fun, v=rbind(N=-1, ct=1))
  initializer <- function(params, t0, ...) {
    c(N=N_0,ct=12)
  }
  pomp(data = data, times = "time", t0 = t0, params = c(mu=mu),
    rprocess = rprocess, initializer = initializer,
    zeronames="ct", paramnames=c("mu"), statenames=c("N","ct"),
    measurement.model = mmodel)
}

simulate(create_example(times = 1), as.data.frame=TRUE)
simulate(create_example(times = c(1,2)), as.data.frame=TRUE)

dev.off()

