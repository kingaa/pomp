options(digits=3)
png(filename="gillespie-%02d.png",res=100)

library(pomp)
library(magrittr)
library(reshape2)
library(ggplot2)

set.seed(754646834L)

c(gamma=24,mu=1/70,iota=0.1,
  beta1=330,beta2=410,beta3=490,
  beta.sd=0,
  rho=0.1,
  S_0=0.07,I_0=1e-4,R_0=0.93,
  pop=1000000
) -> params

cbind(
  birth=c(S=1,I=0,R=0,N=1,cases=0),
  sdeath=c(-1,0,0,-1,0),
  infection=c(-1,1,0,0,0),
  ideath=c(0,-1,0,-1,0),
  recovery=c(0,-1,1,0,1),
  rdeath=c(0,0,-1,-1,0)
) -> Vmatrix

rate.fun <- function(j, N, S, I, R,
  iota, mu, gamma, beta1, beta2, beta3,
  seas_1, seas_2, seas_3, ...) {
  switch(
    j,
    mu*N,            # birth
    mu*S,            # susceptible death
    {                # infection
      beta <- beta1*seas_1+beta2*seas_2+beta3*seas_3
      (beta*I+iota)*S/N
    },
    mu*I,            # infected death
    gamma*I,         # recovery
    mu*R,            # recovered death
    stop("unrecognized event ",j)
  ) -> r
  r
}

simulate(
  params=params,
  seed=806867104L,
  times=seq(from=0,to=2,by=1/52),
  t0=0,
  rprocess=gillespie.sim(
    rate.fun=rate.fun,
    v=Vmatrix,
    hmax=1/52/10
  ),
  zeronames=c("cases"),
  covar=covariate_table(
    t=seq(0,2,by=1/52/10),
    seas=periodic.bspline.basis(
      t,degree=3,period=1,nbasis=3),
    times="t"
  ),
  rmeasure=function(cases,rho,...)
    c(reports=rbinom(n=1,size=round(cases),prob=rho)),
  dmeasure=function(reports,cases,rho,...,log)
    dbinom(x=reports,size=round(cases),prob=rho,log=log),
  rinit=function (S_0, I_0, R_0, pop, ...) {
    s <- S_0+I_0+R_0
    c(N=pop,cases=0,round(pop/s*c(S=S_0,I=I_0,R=R_0)))
  }
) -> gsir

gsir %>%
  simulate(zeronames=NULL) %>%
  plot(main="Gillespie SIR, no zeroing")

gsir %>%
  simulate(
    rprocess=gillespie.sim(
      rate.fun=Csnippet("
        double beta;
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
        beta = dot_product(nbasis,&seas_1,&beta1);
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
      ),
      v=Vmatrix,
      hmax=1/52/10
    ),
    nbasis=3L,
    seed=806867104L,
    paramnames=c("gamma","mu","iota","beta1","beta2","beta3","beta.sd",
      "pop","rho"),
    statenames=c("S","I","R","N","cases")
  ) -> gsir1

gsir1 %>%
  plot(main="Gillespie SIR, with zeroing")

pompExample(sir2)
gsir2 <- simulate(sir2,params=coef(gsir),
  times=time(gsir),t0=timezero(gsir),seed=806867104L)

list(R=as.data.frame(gsir),
  Csnippet=as.data.frame(gsir1),
  sir2=as.data.frame(gsir2)) %>%
  melt(id="time") %>%
  subset(variable=="reports") %>%
  ggplot(aes(x=time,y=value,color=L1))+
  labs(color="",y="reports",title="comparison of implementations")+
  geom_line()+
  theme_bw()+theme(legend.position=c(0.2,0.8))

try(gillespie.sim(rate.fun=rate.fun,v=as.numeric(Vmatrix)))
w <- Vmatrix
colnames(w) <- c(letters[1:5],"a")
try(gillespie.sim(rate.fun=rate.fun,v=w))
w <- Vmatrix
rownames(w) <- c(letters[1:4],"a")
try(gillespie.sim(rate.fun=rate.fun,v=w))
try(gillespie.hl.sim(a="bob"))
try(gillespie.hl.sim(a=list(c("bob","mary"),2)))
try(gillespie.hl.sim(a=list(1,2)))
try(gillespie.hl.sim(.pre=3,a=list("bob",c(a=2,b=1))))
try(gillespie.hl.sim(a=list("bob",c(a="h",b="k"))))
try(gillespie.hl.sim(a=list("bob",c(2,1))))
try(gillespie.hl.sim(a=list("bob",list(a=2,b=1))))
try(gillespie.hl.sim(a=list("bob",c(a=2,1))))
try(gillespie.hl.sim(a=list("bob",c(a=2,b=1)),a=list("mary",c(a=1,a=2))))
try(gillespie.hl.sim(a=list("bob",c(a=2,b=1)),a=list("mary",c(a=1,c=2))))
gillespie.hl.sim(A=list("bob",c(a=2,b=1)),B=list("mary",c(a=1,c=2))) -> f
f@v

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

simulate(gsir,rprocess=gillespie.sim(rate.fun=rate.fun.bad,v=Vmatrix)) %>%
  plot(main="freeze at time 1")

rate.fun.bad <- function(j, x, t, params, covars, ...) {
  if (t>0.1) {
    -rate.fun(j,x,t,params,covars,...)
  } else {
    rate.fun(j,x,t,params,covars,...)
  }
}

try(simulate(gsir,rprocess=gillespie.sim(rate.fun=rate.fun.bad,v=Vmatrix)))

rate.fun.bad <- function(j, x, t, params, covars, ...) -1

try(pomp(gsir,rprocess=gillespie.sim(rate.fun=rate.fun.bad,v=Vmatrix)) %>% simulate())

rate.fun.bad <- function(j, x, t, params, covars, ...) c(1,1)

try(pomp(gsir,rprocess=gillespie.sim(rate.fun=rate.fun.bad,v=Vmatrix)) %>% simulate())

try(pomp(gsir,rprocess=gillespie.sim(rate.fun=3,v=Vmatrix)))

create_example <- function(times = c(1,2), t0 = 0, mu = 0.001, N_0 = 1) {
  rate.fun <- function(j, mu, N, ...) {
    switch(j, mu*N, stop("unrecognized event ",j))
  }
  rprocess <- gillespie.sim(rate.fun = rate.fun, v=rbind(N=-1, ct=1))
  initializer <- function(t0, ...) {
    c(N=N_0,ct=12)
  }
  simulate(times = times, t0 = t0, params = c(mu=mu),
    rprocess = rprocess, rinit = initializer,
    zeronames="ct", format="data.frame")
}

create_example(times = 1)
create_example(times = c(1,2))

dev.off()
