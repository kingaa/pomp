library(pomp)
library(ggplot2)
library(plyr)
library(reshape2)
library(magrittr)

set.seed(49596868L)

png(filename="dp-%02d.png",res=100)

create_example <- function(times = 1, t0 = 0, mu = 0.001, N_0 = 1,
                           simulator = c("gillespie","euler","onestep")) {
    data <- data.frame(time = times, reports = NA)
    d <- cbind(death = c(1,0))
    v <- cbind(death = c(-1,1))
    e <- c(0.03,0)
    f <- function(j, x, t, params, ...){
        params["mu"] * x[1]
    }
    simulator <- match.arg(simulator)
    switch(
        simulator,
        gillespie=gillespie.sim(rate.fun = f, v = v, d = d),
        euler=euler.sim(
            Csnippet("double x = rbinom(N,1-exp(-mu*dt)); N -= x; ct += x;"),
            delta.t=0.1
        ),
        onestep=onestep.sim(
            Csnippet("double x = rbinom(N,1-exp(-mu*dt)); N -= x; ct += x;")
        )
    ) -> rprocess
    initializer <- function(params, t0, ...) {
        c(N=N_0,ct=12)
    }
    pomp(data = data, times = "time", t0 = t0, params = c(mu=mu),
         rprocess = rprocess, initializer = initializer, zeronames="ct",
         paramnames=c("mu"), statenames=c("N","ct"))
}

create_example(simulator="gillespie",times=c(0,1,10,100,1000)) %>%
    simulate(as.data.frame=TRUE, states=TRUE, nsim = 1000) %>%
    count(~time+N) %>%
    ddply(~time,mutate,freq=freq/sum(freq)) %>%
    dcast(time~N,value.var="freq")
create_example(times=seq(0,5,by=0.2),mu=0.01,N_0=100) %>%
    simulate(nsim=100,as.data.frame=TRUE, states=TRUE) -> sims
sims %>%
    subset(sim<=4) %>%
    melt(id=c("time","sim")) %>%
    ggplot(aes(x=time,y=value,group=interaction(sim,variable)))+
    geom_step()+
    facet_grid(variable~sim,scales="free_y")+
    labs(title="death process, Gillespie",subtitle=expression(mu==0.01))
stopifnot(sims %>% ddply(~sim,mutate,s=cumsum(ct),Nn=(N+s)==100) %>%
          subset(!Nn) %>% nrow() %>% equals(0))

create_example(simulator="onestep",times=c(0,1,10,100,1000)) %>%
    simulate(as.data.frame=TRUE, states=TRUE, nsim = 1000) %>%
    count(~time+N) %>%
    ddply(~time,mutate,freq=freq/sum(freq)) %>%
    dcast(time~N,value.var="freq")
create_example(simulator="onestep",
               times=seq(0,5,by=0.2),mu=0.01,N_0=100) %>%
    simulate(nsim=100,as.data.frame=TRUE, states=TRUE) -> sims
sims %>%
    subset(sim<=4) %>%
    melt(id=c("time","sim")) %>%
    ggplot(aes(x=time,y=value,group=interaction(sim,variable)))+
    geom_step()+
    facet_grid(variable~sim,scales="free_y")+
    labs(title="death process, onestep",subtitle=expression(mu==0.01))
stopifnot(sims %>% ddply(~sim,mutate,s=cumsum(ct),Nn=(N+s)==100) %>%
          subset(!Nn) %>% nrow() %>% equals(0))

create_example(simulator="euler",times=c(0,1,10,100,1000)) %>%
    simulate(as.data.frame=TRUE, states=TRUE, nsim = 1000) %>%
    count(~time+N) %>%
    ddply(~time,mutate,freq=freq/sum(freq)) %>%
    dcast(time~N,value.var="freq")
create_example(simulator="euler",
               times=seq(0,5,by=0.2),mu=0.01,N_0=100) %>%
    simulate(nsim=100,as.data.frame=TRUE, states=TRUE) -> sims
sims %>%
    subset(sim<=4) %>%
    melt(id=c("time","sim")) %>%
    ggplot(aes(x=time,y=value,group=interaction(sim,variable)))+
    geom_step()+
    facet_grid(variable~sim,scales="free_y")+
    labs(title="death process, Euler",subtitle=expression(mu==0.01))
stopifnot(sims %>% ddply(~sim,mutate,s=cumsum(ct),Nn=(N+s)==100) %>%
          subset(!Nn) %>% nrow() %>% equals(0))

create_example(mu=1) %>%
  simulate(as.data.frame=TRUE, states=TRUE, times = c(1), nsim = 1000, seed=1066) %>%
  count(~N)
create_example(mu=1) %>%
  simulate(as.data.frame=TRUE, states=TRUE, times = c(0,1), nsim = 1000, seed=1066) %>%
  subset(time>0) %>%
  count(~N)
create_example() %>%
  simulate(as.data.frame=TRUE, states=TRUE, times = c(1e4), nsim = 10000, seed=1066) %>%
  count(~N)

create_example(N_0=1000,mu=0.02,simulator="gillespie",times=c(0,1,10,100,1000)) %>%
  simulate(as.data.frame=TRUE, states=TRUE, nsim = 1000, seed=374244) %>%
  ggplot(aes(x=N,group=time))+
  geom_histogram(aes(y=..density..),binwidth=10)+
  facet_grid(time~.)+
  labs(title="death process, Gillespie",subtitle=expression(mu==0.02))

create_example(N_0=1000,mu=0.02,simulator="onestep",times=c(0,1,10,100,1000)) %>%
  simulate(as.data.frame=TRUE, states=TRUE, nsim = 1000, seed=374244) %>%
  ggplot(aes(x=N,group=time))+
  geom_histogram(aes(y=..density..),binwidth=10)+
  facet_grid(time~.)+
  labs(title="death process, onestep",subtitle=expression(mu==0.02))

create_example(N_0=1000,mu=0.02,simulator="euler",times=c(0,1,10,100,1000)) %>%
  simulate(as.data.frame=TRUE, states=TRUE, nsim = 1000, seed=374244) %>%
  ggplot(aes(x=N,group=time))+
  geom_histogram(aes(y=..density..),binwidth=10)+
  facet_grid(time~.)+
  labs(title="death process, Euler",subtitle=expression(mu==0.02))

dev.off()
