library(pomp)
library(foreach)
library(doParallel)
library(ggplot2)
library(plyr)
library(reshape2)
library(magrittr)
theme_set(theme_bw())

set.seed(1666441244L,kind="L'Ecuyer-CMRG")

png(filename="parallel-%02d.png",res=100)

registerDoParallel(cores=2)

foreach (i=1:2,.options.multicore=list(set.seed=TRUE)) %dopar% {

    data.frame(time=seq(1,10),P=NA) %>%
        pomp(
            times='time',t0=0,
            rprocess=euler.sim(
                Csnippet("
              double dW = rnorm(0,sqrt(dt));
              N += r*N*(1-N/K)*dt+sigma*N*dW;"),
              delta.t=0.05),
            rmeasure=Csnippet("P = rpois(N);"),
            dmeasure=Csnippet("lik = dpois(P,N,give_log);"),
            skeleton=vectorfield(Csnippet("DN = r*N*(1-N/K);")),
            initializer=Csnippet("N = exp(rnorm(log(N_0),1));"),
            paramnames=c("sigma","r","K","N_0"),
            statenames="N",
            params=c(r=1,K=200,sigma=0.5,N_0=20)) %>%
        simulate(seed=1323046840L) -> mod

    traj <- trajectory(mod,as.data.frame=TRUE)
    sim <- simulate(mod,nsim=5,states=TRUE,as.data.frame=TRUE)
    list(trajectory=traj,simulations=sim)
} %>% melt(measure.vars="N") %>%
    dcast(L2+variable+time~L1+sim+traj) %>%
    subset(select=-variable) %>%
    melt(id=c("L2","time")) %>% 
    ggplot(aes(x=time,y=value,group=interaction(variable,L2),
               color=L2))+
    geom_line()+
    expand_limits(y=c(0,NA))+
    labs(color="",y="N")+
    theme(legend.position=c(0.2,0.9))

dev.off()
