library(pomp)
library(magrittr)
library(reshape2)
library(ggplot2)

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

cbind(
    birth=c(1,0,0,1,0),
    sdeath=c(-1,0,0,-1,0),
    infection=c(-1,1,0,0,0),
    ideath=c(0,-1,0,-1,0),
    recovery=c(0,-1,1,0,1),
    rdeath=c(0,0,-1,-1,0)
) -> Vmatrix

cbind(
    birth=c(0,0,0,1,0),
    sdeath=c(1,0,0,0,0),
    infection=c(1,1,0,1,0),
    ideath=c(0,1,0,0,0),
    recovery=c(0,1,0,0,0),
    rdeath=c(0,0,1,0,0)
) -> Dmatrix

data.frame(
    time=seq(from=0,to=2,by=1/52),
    reports=NA
) %>% 
    pomp(
        times="time",
        t0=0,
        rprocess=gillespie.sim(rate.fun=rate.fun,v=Vmatrix,d=Dmatrix),
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

pompExample(gillespie.sir)
gsir2 <- simulate(gillespie.sir,params=coef(gsir),
                  times=time(gsir),t0=timezero(gsir),seed=806867104L)

list(R=as.data.frame(gsir),
     C=as.data.frame(gsir2)) %>%
    melt(id="time") %>%
    subset(variable=="reports") %>%
    ggplot(aes(x=time,y=value,color=L1))+
    labs(color="",y="reports")+
    geom_line()+
    theme_bw()+theme(legend.position=c(0.2,0.8)) -> pl1

gsir %>%
    pomp(
        rprocess=kleap.sim(
            rate.fun=rate.fun,
            e=rep(0.1,5),
            v=Vmatrix,
            d=Dmatrix),
    ) -> gsir2

simulate(gsir,seed=95424809L,as.data.frame=TRUE) -> exact
simulate(gsir2,seed=95424809L,as.data.frame=TRUE) -> kleap

list(exact=exact,`K leap`=kleap) %>%
    melt(id=c("time","sim")) %>%
    subset(variable=="cases") %>%
    ggplot(aes(x=time,y=value,group=interaction(sim,L1),color=L1))+
    labs(color="",y="cases")+
    geom_line()+
    theme_bw()+theme(legend.position=c(0.2,0.8)) -> pl2

rate.fun.bad <- function(j, x, t, params, covars, ...) {
    -rate.fun(j,x,t,params,covars,...)
}

try(pomp(gsir,rprocess=gillespie.sim(rate.fun=rate.fun.bad,v=Vmatrix,d=Dmatrix)) %>% simulate())
try(pomp(gsir,rprocess=kleap.sim(rate.fun=rate.fun.bad,e=rep(0.2,5),v=Vmatrix,d=Dmatrix)) %>% simulate())

rate.fun.bad <- function(j, x, t, params, covars, ...) {
    0
}

try(pomp(gsir,rprocess=gillespie.sim(rate.fun=rate.fun.bad,v=Vmatrix,d=Dmatrix)) %>% simulate())
try(pomp(gsir,rprocess=kleap.sim(rate.fun=rate.fun.bad,e=rep(0.2,5),v=Vmatrix,d=Dmatrix)) %>% simulate())

ggsave(plot=pl1,filename="gillespie-01.png",dpi=100,height=4,width=4)
ggsave(plot=pl2,filename="gillespie-02.png",dpi=100,height=4,width=4)
