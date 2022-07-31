params <-
list(prefix = "getting_started")

## ----prelims,echo=FALSE,cache=FALSE-------------------------------------------
options(
  keep.source=TRUE,
  stringsAsFactors=FALSE,
  encoding="UTF-8",
  scipen=5,
  dplyr.summarise.inform=FALSE,
  pomp_archive_dir="results/getting_started"
)
set.seed(594709947L)
library(tidyverse)
theme_set(theme_bw())
bigtick <- Sys.time()


## ----parallel_setup,include=FALSE,purl=TRUE,cache=FALSE-----------------------
if (file.exists("CLUSTER.R")) {
  source("CLUSTER.R")
} else {
  library(doParallel)
  registerDoParallel()
}
library(doRNG)
registerDoRNG(348885445L)




## ----sim1,warning=TRUE--------------------------------------------------------
library(pomp)

simulate(t0=0, times=1:20,
  params=c(r=1.2,K=200,sigma=0.1,N_0=50),
  rinit=function (N_0, ...) {
    c(N=N_0)
  },
  rprocess=discrete_time(
    function (N, r, K, sigma, ...) {
      eps <- rnorm(n=1,mean=0,sd=sigma)
      c(N=r*N*exp(1-N/K+eps))
    },
    delta.t=1
  )
) -> sim1


## ----sim1_print---------------------------------------------------------------
sim1




## ----sim1_plot----------------------------------------------------------------
plot(sim1)


## ----sim1_print2--------------------------------------------------------------
as(sim1,"data.frame")


## ----sim1_plot2---------------------------------------------------------------
  ggplot(data=as.data.frame(sim1),aes(x=time,y=N))+
  geom_line()


## ----sim2---------------------------------------------------------------------
simulate(t0=0, times=1:20,
  params=c(r=1.2,K=200,sigma=0.1,N_0=50,b=0.05),
  rinit=function (N_0, ...) {
    c(N=N_0)
  },
  rprocess=discrete_time(
    function (N, r, K, sigma, ...) {
      eps <- rnorm(n=1,mean=0,sd=sigma)
      c(N=r*N*exp(1-N/K+eps))
    },
    delta.t=1
  ),
  rmeasure=function (N, b, ...) {
    c(Y=rpois(n=1,lambda=b*N))
  }
) -> sim2


## ----sim2_alt-----------------------------------------------------------------
simulate(
  sim1,
  params=c(r=1.2,K=200,sigma=0.1,N_0=50,b=0.05),
  rmeasure=function (N, b, ...) {
    c(Y=rpois(n=1,lambda=b*N))
  }
) -> sim2



## ----sim2_print_plot----------------------------------------------------------
as(sim2,"data.frame")
plot(sim2)
ggplot(data=gather(
  as(sim2,"data.frame"),
  variable,value,-time),
  aes(x=time,y=value,color=variable))+
  geom_line()


## ----sim2_sim1----------------------------------------------------------------
simulate(sim2,nsim=20) -> sims

ggplot(data=gather(
  as.data.frame(sims),
  variable,value,Y,N
),
  aes(x=time,y=value,color=variable,
    group=interaction(.id,variable)))+
  geom_line()+
  facet_grid(variable~.,scales="free_y")+
  labs(y="",color="")


## ----sim2_sim2----------------------------------------------------------------
p <- parmat(coef(sim2),3)
p["sigma",] <- c(0.05,0.25,1)
colnames(p) <- LETTERS[1:3]

simulate(sim2,params=p,format="data.frame") -> sims
sims <- gather(sims,variable,value,Y,N)
ggplot(data=sims,aes(x=time,y=value,color=variable,
  group=interaction(.id,variable)))+
  geom_line()+
  scale_y_log10()+
  expand_limits(y=1)+
  facet_grid(variable~.id,scales="free_y")+
  labs(y="",color="")


## ----sim2_sim3,message=F------------------------------------------------------
simulate(sim2,params=p,
  times=seq(0,3),
  nsim=500,format="data.frame") -> sims
ggplot(data=separate(sims,.id,c("parset","rep")),
  aes(x=N,fill=parset,group=parset,color=parset))+
  geom_density(alpha=0.5)+
  # geom_histogram(aes(y=..density..),position="dodge")+
  facet_grid(time~.,labeller=label_both,scales="free_y")+
  lims(x=c(NA,1000))


## ----parus-plot---------------------------------------------------------------
parus |>
  ggplot(aes(x=year,y=pop))+
  geom_line()+geom_point()+
  expand_limits(y=0)
parus


## ----rick---------------------------------------------------------------------
parus |>
  pomp(
    times="year", t0=1960,
    rinit=function (N_0, ...) {
      c(N=N_0)
    },
    rprocess=discrete_time(
      function (N, r, K, sigma, ...) {
        eps <- rnorm(n=1,mean=0,sd=sigma)
        c(N=r*N*exp(1-N/K+eps))
      },
      delta.t=1
    ),
    rmeasure=function (N, b, ...) {
      c(pop=rpois(n=1,lambda=b*N))
    }
  ) -> rick


## ----logistic-step-fun--------------------------------------------------------
vpstep <- function (N, r, K, sigma, delta.t, ...) {
  dW <- rnorm(n=1,mean=0,sd=sqrt(delta.t))
  c(N = N + r*N*(1-N/K)*delta.t + sigma*N*dW)
}



## ----vp-----------------------------------------------------------------------
rick |> pomp(rprocess=euler(vpstep,delta.t=1/365)) -> vp


## ----vp_sim1------------------------------------------------------------------
vp |>
  simulate(
    params=c(r=0.5,K=2000,sigma=0.1,b=0.1,N_0=2000),
    format="data.frame", include.data=TRUE, nsim=5) |>
  mutate(ds=case_when(.id=="data"~"data",TRUE~"simulation")) |>
  ggplot(aes(x=year,y=pop,group=.id,color=ds))+
  geom_line()+
  labs(color="")


## ----two_snips----------------------------------------------------------------
Csnippet("
  pop = rpois(b*N);  
  ") -> rmeas

Csnippet("
  N = N_0;
  ") -> rinit


## ----more_snips---------------------------------------------------------------
Csnippet("
  double eps = rnorm(0,sigma);
  N = r*N*exp(1-N/K+eps);
") -> rickstepC

Csnippet("
  double dW = rnorm(0,sqrt(dt));
  N += r*N*(1-N/K)*dt+sigma*N*dW;
") -> vpstepC


## ----rickC--------------------------------------------------------------------
parus |>
  pomp(
    times="year", t0=1960,
    rinit=rinit,
    rmeasure=rmeas,
    rprocess=discrete_time(rickstepC,delta.t=1),
    statenames="N",
    paramnames=c("r","K","sigma","b","N_0")
  ) -> rickC


## ----vpC----------------------------------------------------------------------
parus |>
  pomp(
    times="year", t0=1960,
    rinit=rinit,
    rmeasure=rmeas,
    rprocess=euler(vpstepC,delta.t=1/365),
    statenames="N",
    paramnames=c("r","K","sigma","b","N_0")
  ) -> vpC


## ----logistic-dmeasure--------------------------------------------------------
Csnippet("
  lik = dpois(pop,b*N,give_log);
") -> dmeas


## ----pfilter1-----------------------------------------------------------------
rickC |>
  pfilter(Np=1000,
    params=c(r=1.2,K=2000,sigma=0.3,N_0=1600,b=0.1),
    dmeasure=dmeas,
    paramnames="b",statenames="N") -> pfrick


## ----pfrick_print-------------------------------------------------------------
pfrick


## ----pfrick_loglik------------------------------------------------------------
logLik(pfrick)


## ----pfrick_plot--------------------------------------------------------------
plot(pfrick)
as(pfrick,"data.frame")


## ----rick_loglik--------------------------------------------------------------
replicate(10, pfrick |> pfilter() |> logLik()) -> lls
lls
logmeanexp(lls,se=TRUE) -> ll_rick1
ll_rick1



## ----vp_loglik_eval,include=FALSE,purl=TRUE-----------------------------------
bake("vp_loglik.rds",{
  replicate(10,
    vpC |>
      pfilter(Np=1000,dmeasure=dmeas,
        params=c(r=0.5,K=2000,sigma=0.1,b=0.1,N_0=2000),
        paramnames="b",statenames="N") |>
      logLik()) |>
    logmeanexp(se=TRUE) -> ll_vp1
  ll_vp1
}) -> ll_vp1
ll_vp1


## ----rick_skel----------------------------------------------------------------
rickC |>
  pomp(
    skeleton=map(
      Csnippet("DN = r*N*exp(1-N/K);"),
      delta.t=1
    ),
    paramnames=c("r","K"), statenames="N"
  ) -> rickC


## ----vp_skel------------------------------------------------------------------
vpC |>
  pomp(
    skeleton=vectorfield(Csnippet("DN = r*N*(1-N/K);")),
    paramnames=c("r","K"), statenames="N"
  ) -> vpC


## ----vp_traj1-----------------------------------------------------------------
p <- parmat(c(r=0.5,K=2000,sigma=0.1,b=0.1,N_0=2000),10)
p["N_0",] <- seq(10,3000,length=10)
vpC |>
  trajectory(params=p,format="data.frame") |>
  ggplot(mapping=aes(x=year,y=N,color=.id,group=.id))+
  guides(color="none")+
  geom_line()+
  theme_bw()


## ----vp_traj_objfun-----------------------------------------------------------
vpC |>
  traj_objfun(
    est=c("K","N_0"),
    params=c(r=0.5,K=2000,sigma=0.1,b=0.1,N_0=2000),
    dmeasure=dmeas, statenames="N", paramnames="b"
  ) -> ofun


## ----ofun_print---------------------------------------------------------------
ofun


## ----ofun_eval----------------------------------------------------------------
ofun(c(1000,3000))


## ----traj_match1--------------------------------------------------------------
library(subplex)

subplex(c(2000,1500),fn=ofun) -> fit
fit 


## ----traj_match_print---------------------------------------------------------
ofun(fit$par)
coef(ofun)
logLik(ofun)


## ----traj_match_traj----------------------------------------------------------
ofun |>
  trajectory(format="data.frame") |>
  ggplot(mapping=aes(x=year,y=N,color=.id,group=.id))+
  guides(color="none")+
  geom_line()+
  theme_bw()


## ----traj_match_plot----------------------------------------------------------
ofun |>
  trajectory(format="data.frame") |>
  mutate(
    pop=coef(ofun,"b")*N,
    .id="prediction"
  ) |>
  select(-N) |>
  rbind(
    parus |>
      mutate(
        .id="data",
        pop=as.double(pop)
      )
  ) |>
  spread(.id,pop) |>
  ggplot(aes(x=year))+
  geom_line(aes(y=prediction))+
  geom_point(aes(y=data))+
  expand_limits(y=0)+
  labs(y="pop")


## ----vp_partrans--------------------------------------------------------------
vpr_partrans <- parameter_trans(log=c("r","K","sigma","b","N_0"))


## ----ofun2--------------------------------------------------------------------
vpC |>
  traj_objfun(
    est=c("b","K","N_0"),
    params=c(r=0.5,K=2000,sigma=0.1,b=0.1,N_0=2000),
    partrans=parameter_trans(log=c("K","N_0","b")),
    dmeasure=dmeas, statenames="N", paramnames=c("b","K","N_0")
  ) -> ofun2

subplex(log(c(0.1,2000,1500)),fn=ofun2) -> fit
ofun2(fit$par)
coef(ofun2)


## ----mif_guesses,fig.width=4,fig.height=4-------------------------------------
sobol_design(
  lower=c(r=0,K=100,sigma=0,N_0=150,b=1),
  upper=c(r=5,K=600,sigma=2,N_0=150,b=1),
  nseq=100
) -> guesses
guesses
plot(guesses,pch=16)



## ----vp_mif1_eval,include=FALSE,purl=TRUE-------------------------------------
bake("vp_mif1.rds",{
  vpC |>
    mif2(
      params=guesses[1,],
      Np=1000,
      Nmif=20,
      dmeasure=dmeas,
      partrans=parameter_trans(log=c("r","K","sigma","N_0")),
      rw.sd=rw.sd(r=0.02,K=0.02,sigma=0.02,N_0=ivp(0.02)),
      cooling.fraction.50=0.5,
      paramnames=c("r","K","sigma","N_0","b"),
      statenames=c("N")
    ) -> mf1
}) -> mf1


## ----mf1_display--------------------------------------------------------------
mf1


## ----mif_plot1----------------------------------------------------------------
plot(mf1)


## ----mf_pfilter1--------------------------------------------------------------
replicate(5, mf1 |> pfilter() |> logLik()) |> logmeanexp(se=TRUE)



## ----parus_mif1_eval,include=FALSE,eval=TRUE,purl=TRUE------------------------
bake("parus_mif1.rds",{
  library(foreach)
  
  foreach (guess=iter(guesses,"row"),
    .combine=c, .packages=c("pomp"),
    .errorhandling="remove", .inorder=FALSE) %dopar% {
      
      mf1 |> mif2(params=guess)
      
    } -> mifs
}) -> mifs


## ----mif_plot2----------------------------------------------------------------
mifs |>
  traces() |>
  melt() |>
  filter(variable!="b") |>
  ggplot(aes(x=iteration,y=value,group=L1,color=L1))+
  geom_line()+
  facet_wrap(~variable,scales="free_y")+
  guides(color="none")



## ----mif_plot3----------------------------------------------------------------
mifs |> 
  as("data.frame") |> 
  gather(variable,value,-year,-.id) |>
  ggplot(aes(x=year,y=value,group=.id,color=.id))+
  geom_line()+
  facet_wrap(~variable,scales="free_y",ncol=1)+
  guides(color="none")



## ----parus_pf1_eval,include=FALSE,eval=TRUE,purl=TRUE-------------------------
bake(file="parus_pf1.rds",{
  foreach (mf=mifs,
    .combine=rbind, .packages=c("pomp"), 
    .errorhandling="remove", .inorder=FALSE) %dopar% {
      
      replicate(5, 
        mf |> pfilter() |> logLik()
      ) |>
        logmeanexp(se=TRUE) -> ll
      
      data.frame(as.list(coef(mf)),loglik=ll[1],loglik.se=ll[2])
      
    } -> estimates
}) -> estimates


## ----mif_plot4,warning=FALSE,fig.width=4,fig.height=4-------------------------
estimates |>
  full_join(guesses) |>
  filter(is.na(loglik) | loglik>max(loglik,na.rm=TRUE)-30) |>
  do({
    pairs(~loglik+r+K+sigma+N_0,data=.,
      pch=16,
      col=if_else(is.na(.$loglik),"#99999955","#ff0000ff"))
    .
  }) |> 
  invisible()





## ----parus_mif2_eval,include=FALSE,eval=TRUE,purl=TRUE------------------------
bake(file="parus_mif2.rds",{
  estimates |>
    filter(!is.na(loglik)) |>
    filter(loglik > max(loglik)-30) |>
    select(-loglik,-loglik.se) -> starts
  
  foreach (start=iter(starts,"row"),
    .combine=rbind, .packages=c("pomp"), 
    .errorhandling="remove", .inorder=FALSE) %dopar% {
      
      mf1 |> 
        mif2(params=start) |>
        mif2() -> mf
      
      replicate(5, 
        mf |> pfilter() |> logLik()
      ) |>
        logmeanexp(se=TRUE) -> ll
      
      data.frame(as.list(coef(mf)),loglik=ll[1],loglik.se=ll[2])
      
    } -> ests1
}) -> ests1


## ----db1----------------------------------------------------------------------
estimates |>
  rbind(ests1) -> estimates

estimates |>
  arrange(-loglik)


## ----mif2_plot1,fig.width=4,fig.height=4--------------------------------------
estimates |>
  filter(loglik>max(loglik,na.rm=TRUE)-4) |>
  do({
    pairs(~loglik+r+K+sigma+N_0,data=.,
      pch=16,
      col=if_else(is.na(.$loglik),"#99999955","#ff0000ff"))
    .
  }) |> 
  invisible()


## ----parus_pd,fig.width=4,fig.height=4----------------------------------------
estimates |>
  filter(loglik>max(loglik)-10) |>
  select(r,K,sigma,N_0,b) |>
  apply(2,range) -> ranges
ranges

profile_design(
  r=10^seq(
    from=log10(ranges[1,1]),
    to=log10(ranges[2,1]),
    length=20
  ),
  lower=ranges[1,-1],
  upper=ranges[2,-1],
  nprof=50
) -> starts

dim(starts)

pairs(~r+K+sigma+N_0+b,data=starts)


## ----parus_profile_eval,include=FALSE,purl=TRUE,eval=TRUE---------------------
bake(file="parus_profile.rds",{
  foreach (start=iter(starts,"row"),
    .combine=rbind, .packages=c("pomp"),
    .errorhandling="remove", .inorder=FALSE) %dopar% {
    
    mf1 |>
      mif2(
        params=start,
        partrans=parameter_trans(log=c("K","sigma","N_0")),
        rw.sd=rw.sd(K=0.02,sigma=0.02,N_0=ivp(0.02)),
        paramnames=c("K","sigma","N_0","b")
      ) |>
      mif2() -> mf
    
    replicate(5, 
      mf |> pfilter() |> logLik()
    ) |>
      logmeanexp(se=TRUE) -> ll
    
    data.frame(as.list(coef(mf)),loglik=ll[1],loglik.se=ll[2])
  } -> r_prof
}) -> r_prof


## ----db2----------------------------------------------------------------------
estimates |>
  rbind(r_prof) -> estimates


## ----prof_plot1---------------------------------------------------------------
r_prof |>
  group_by(r) |>
  filter(rank(-loglik)<=2) |>
  ungroup() |>
  ggplot(aes(x=r,y=loglik,
    ymin=loglik-2*loglik.se,ymax=loglik+2*loglik.se))+
  geom_point()+
  geom_errorbar()+
  geom_smooth(method="loess",span=0.2)+
  scale_x_log10()


## ----r_prof_lrt,echo=FALSE----------------------------------------------------
r_prof |>
  group_by(r) |>
  filter(rank(-loglik)<=2) |>
  ungroup() |>
  do({
    loess(loglik~log(r),data=.,span=0.25) -> fit
    rr <- range(.$r)
    r <- exp(seq(log(rr[1]),log(rr[2]),length=200))
    data.frame(r=r,loglik=predict(fit,newdata=log(r)))
  }) -> prof

crit <- 0.5*qchisq(df=1,p=0.95)
cutoff <- max(prof$loglik)-crit
ci <- range(prof$r[prof$loglik>cutoff])


## ----prof_plot2---------------------------------------------------------------
r_prof |>
  group_by(r) |>
  filter(rank(-loglik)<=2) |>
  ungroup() |>
  ggplot(aes(x=r,y=sigma))+
  geom_point()+
  geom_smooth(method="loess",span=0.2)+
  scale_x_log10()+
  labs(y=expression(sigma))


## ----parus-pmcmc-starts-------------------------------------------------------
r_prof |>
  group_by(r) |>
  filter(loglik==max(loglik)) |>
  ungroup() |>
  filter(
    r > 0.75,
    r < 6,
    sigma < 2,
    K > 100,
    K < 600,
    N_0 > 100,
    N_0 < 600
  ) |>
  select(-loglik,-loglik.se) -> starts
starts




## ----parus-pmcmc-eval,eval=TRUE,purl=TRUE,include=FALSE-----------------------
bake(file="parus-pmcmc.rds",dependson=starts,{
  foreach (start=iter(starts,"row"),.combine=c,
    .errorhandling="remove",.inorder=FALSE) %dopar% 
    {
      library(pomp)
      mf1 |>
        pmcmc(
          Nmcmc=2000,Np=200,params=start,
          dprior=Csnippet("
            lik = dunif(r,0,10,1)+dunif(sigma,0,2,1)+
                  dunif(K,0,600,1)+dunif(N_0,0,600,1);
            lik = (give_log) ? lik : exp(lik);"
          ),
          paramnames=c("K","N_0","r","sigma"),
          proposal=mvn.rw.adaptive(
            rw.sd=c(r=0.02,sigma=0.02,K=50,N_0=50),
            scale.start=100,shape.start=100
          )
        ) -> chain
      chain |> pmcmc(Nmcmc=20000,proposal=mvn.rw(covmat(chain)))
    } -> chains
}) -> chains


## ----pmcmc-diagnostics1-------------------------------------------------------
library(coda)
chains |> traces() -> traces
rejectionRate(traces[,c("r","sigma","K","N_0")])


## ----pmcmc-diagnostics2-------------------------------------------------------
traces %>% autocorr.diag(lags=c(1,5,10,50,100))
traces <- window(traces,thin=100,start=2000)


## ----pmcmc-diagnostics3,fig.dim=c(8,6),out.width="95%"------------------------
traces %>%
  lapply(as.data.frame) %>%
  lapply(rownames_to_column,"iter") %>%
  bind_rows(.id="chain") %>%
  mutate(iter=as.numeric(iter)) %>%
  select(chain,iter,loglik,r,sigma,K,N_0) %>%
  pivot_longer(c(-chain,-iter)) %>%
  ggplot(aes(x=iter,group=chain,color=chain,y=value))+
  guides(color="none")+
  labs(x="iteration",y="")+
  geom_line(alpha=0.3)+geom_smooth(method="loess",se=FALSE)+
  facet_wrap(name~.,scales="free_y",strip.position="left",ncol=2)+
  theme(
    strip.placement="outside",
    strip.background=element_rect(fill=NA,color=NA)
  )

gelman.diag(traces[,c("r","sigma","K","N_0")])


## ----parus_posterior----------------------------------------------------------
traces %>%
  lapply(as.data.frame) %>%
  lapply(rownames_to_column,"iter") %>%
  bind_rows(.id="chain") %>%
  select(chain,iter,loglik,r,sigma,K,N_0) %>%
  pivot_longer(c(-chain,-iter)) %>%
  ggplot(aes(x=value))+
  geom_density()+
  geom_rug()+
  labs(x="")+
  facet_wrap(~name,scales="free",strip.position="bottom")+
  theme(
    strip.placement="outside",
    strip.background=element_rect(fill=NA,color=NA)
  )

traces %>% summary()


## ----mle_sim_plot1------------------------------------------------------------
r_prof |> 
  filter(loglik==max(loglik)) -> mle

mlepomp <- as(mifs[[1]],"pomp")
coef(mlepomp) <- mle

mlepomp |>
  simulate(nsim=8,format="data.frame",include.data=TRUE) |>
  ggplot(mapping=aes(x=year,y=pop,group=.id,alpha=(.id=="data")))+
  scale_alpha_manual(values=c(`TRUE`=1,`FALSE`=0.2),
    labels=c(`FALSE`="simulation",`TRUE`="data"))+
  labs(alpha="")+
  geom_line()+
  theme_bw()


## ----mle_sim_plot2------------------------------------------------------------
mlepomp |>
  simulate(nsim=11,format="data.frame",include.data=TRUE) |>
  ggplot(mapping=aes(x=year,y=pop,group=.id,color=(.id=="data")))+
  scale_color_manual(values=c(`TRUE`="black",`FALSE`="grey50"),
    labels=c(`FALSE`="simulation",`TRUE`="data"))+
  labs(color="")+
  geom_line()+
  facet_wrap(~.id)+
  theme_bw()


## ----probe1-------------------------------------------------------------------
mlepomp |>
  probe(nsim=200,probes=list(
    mean=probe.mean("pop"),
    q=probe.quantile("pop",probs=c(0.05,0.25,0.5,0.75,0.95)),
    probe.acf("pop",lags=c(1,3),type="corr",transform=log)
  )) -> vp_probe

vp_probe


## ----probe1_summary-----------------------------------------------------------
summary(vp_probe)


## ----probe1_plot,fig.width=6.8,fig.height=6.8---------------------------------
plot(vp_probe)


## ----sessinfo-----------------------------------------------------------------
bigtock <- Sys.time()
totalSweaveTime <- bigtock-bigtick
sysi <- Sys.info()
sess <- sessionInfo()
tfile <- file.path("results","getting_started","timing.rda")

if (file.exists(tfile)) {
  load(tfile)
} else {
  save(totalSweaveTime,sysi,sess,file=tfile,compress='xz')
}

print(sysi)
print(sess)
print(totalSweaveTime)

