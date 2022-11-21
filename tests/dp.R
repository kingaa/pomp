options(digits=3)
png(filename="dp-%02d.png",res=100)

library(pomp)
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

set.seed(49596868L)

create_example <- function(times = 1, t0 = 0, mu = 0.001, N_0 = 1,
  simulator = c("gillespie","euler","onestep")) {

  v <- cbind(death = c(-1,1))
  simulator <- match.arg(simulator)
  switch(
    simulator,
    gillespie=gillespie(Csnippet("rate = mu * N;"), v = v),
    euler=euler(
      Csnippet("double x = rbinom(N,1-exp(-mu*dt)); N -= x; ct += x;"),
      delta.t=0.1
    ),
    onestep=onestep(
      Csnippet("double x = rbinom(N,1-exp(-mu*dt)); N -= x; ct += x;")
    )
  ) -> rprocess

  rinit <- Csnippet("N = N_0; ct = 12;")

  pomp(data=NULL, times=times, t0=t0,
    params=c(mu=mu,N_0=N_0),
    rprocess=rprocess, rinit=rinit, accumvars="ct",
    paramnames=c("mu","N_0"), statenames=c("N","ct"))
}

create_example(simulator="gillespie",times=c(0,1,10,100,1000)) %>%
  simulate(format="data.frame", nsim=1000) %>%
  count(time,N) %>%
  group_by(time) %>%
  mutate(freq=n/sum(n)) %>%
  ungroup() %>%
  select(-n) %>%
  pivot_wider(names_from=N,values_from=freq) %>%
  as.data.frame()
create_example(times=seq(0,5,by=0.2),mu=0.01,N_0=100) %>%
  simulate(nsim=100,format="data.frame") -> sims
sims %>%
  filter(.id<=4) %>%
  melt(id=c("time",".id")) %>%
  ggplot(aes(x=time,y=value,group=interaction(.id,variable)))+
  geom_step()+
  facet_grid(variable~.id,scales="free_y")+
  labs(title="death process, Gillespie",subtitle=expression(mu==0.01))
stopifnot(
(sims %>%
   group_by(.id) %>%
   mutate(s=cumsum(ct),Nn=(N+s)==100) %>%
   filter(!Nn) %>%
   nrow()) == 0
)

create_example(simulator="onestep",times=c(0,1,10,100,1000)) %>%
  simulate(format="data.frame", nsim=1000) %>%
  count(time,N) %>%
  group_by(time) %>%
  mutate(freq=n/sum(n)) %>%
  ungroup() %>%
  select(-n) %>%
  pivot_wider(names_from=N,values_from=freq) %>%
  as.data.frame()
create_example(simulator="onestep",
  times=seq(0,5,by=0.2),mu=0.01,N_0=100) %>%
  simulate(nsim=100,format="data.frame") -> sims
sims %>%
  filter(.id<=4) %>%
  melt(id=c("time",".id")) %>%
  ggplot(aes(x=time,y=value,group=interaction(.id,variable)))+
  geom_step()+
  facet_grid(variable~.id,scales="free_y")+
  labs(title="death process, onestep",subtitle=expression(mu==0.01))
stopifnot(
(sims %>%
   group_by(.id) %>%
   mutate(s=cumsum(ct),Nn=(N+s)==100) %>%
   filter(!Nn) %>%
   nrow()) == 0
)

create_example(simulator="euler",times=c(0,1,10,100,1000)) %>%
  simulate(format="data.frame", nsim=1000) %>%
  count(time,N) %>%
  group_by(time) %>%
  mutate(freq=n/sum(n)) %>%
  ungroup() %>%
  select(-n) %>%
  pivot_wider(names_from=N,values_from=freq) %>%
  as.data.frame()
create_example(simulator="euler",
  times=seq(0,5,by=0.2),mu=0.01,N_0=100) %>%
  simulate(nsim=100,format="data.frame") -> sims
sims %>%
  filter(.id<=4) %>%
  melt(id=c("time",".id")) %>%
  ggplot(aes(x=time,y=value,group=interaction(.id,variable)))+
  geom_step()+
  facet_grid(variable~.id,scales="free_y")+
  labs(title="death process, Euler",subtitle=expression(mu==0.01))
stopifnot(
(sims %>%
   group_by(.id) %>%
   mutate(s=cumsum(ct),Nn=(N+s)==100) %>%
   filter(!Nn) %>%
   nrow()) == 0
)

create_example(mu=1) %>%
  simulate(format="data.frame", times=c(1), nsim=1000, seed=1066) %>%
  count(N) %>%
  as.data.frame()
create_example(mu=1) %>%
  simulate(format="data.frame", times=c(0,1), nsim=1000, seed=1066) %>%
  filter(time>0) %>%
  count(N) %>%
  as.data.frame()
create_example() %>%
  simulate(format="data.frame", times=c(1e4), nsim=10000, seed=1066) %>%
  count(N) %>%
  as.data.frame()

create_example(N_0=1000,mu=0.02,simulator="gillespie",
  times=-1/0.02*log(c(1,0.8,0.6,0.4,0.2,0.01))) %>%
  simulate(format="data.frame", nsim=1000, seed=374244) %>%
  ggplot(aes(x=N,group=time))+
  geom_histogram(aes(y=after_stat(density)),binwidth=10)+
  labs(title="death process, Gillespie",subtitle=expression(mu==0.02))+
  facet_grid(time~.,labeller=label_bquote(t==.(time)))+
  theme(strip.text=element_text(size=6))

create_example(N_0=1000,mu=0.02,simulator="onestep",
  times=-1/0.02*log(c(1,0.8,0.6,0.4,0.2,0.01))) %>%
  simulate(format="data.frame", nsim=1000, seed=374244) %>%
  ggplot(aes(x=N,group=time))+
  geom_histogram(aes(y=after_stat(density)),binwidth=10)+
  labs(title="death process, onestep",subtitle=expression(mu==0.02))+
  facet_grid(time~.,labeller=label_bquote(t==.(time)))+
  theme(strip.text=element_text(size=6))

create_example(N_0=1000,mu=0.02,simulator="euler",
  times=-1/0.02*log(c(1,0.8,0.6,0.4,0.2,0.01))) %>%
  simulate(format="data.frame", nsim=1000, seed=374244) %>%
  ggplot(aes(x=N,group=time))+
  geom_histogram(aes(y=after_stat(density)),binwidth=10)+
  labs(title="death process, Euler",subtitle=expression(mu==0.02))+
  facet_grid(time~.,labeller=label_bquote(t==.(time)))+
  theme(strip.text=element_text(size=6))

dev.off()
