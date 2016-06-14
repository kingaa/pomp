## ----prelims,echo=FALSE,cache=FALSE--------------------------------------
library(pomp)
options(
  keep.source=TRUE,
  stringsAsFactors=FALSE,
  encoding="UTF-8",
  scipen=5
)
set.seed(594709947L)
pdf(file="getting_started.pdf")

## ----prelims2,echo=FALSE,cache=FALSE-------------------------------------
library(ggplot2)
library(plyr)
library(reshape2)
library(magrittr)
theme_set(theme_bw())

## ----load-parus-data-----------------------------------------------------
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

## ----parus-plot----------------------------------------------------------
ggplot(data=parus.dat,mapping=aes(x=year,y=P))+
  geom_line()+geom_point()+
  expand_limits(y=0)+
  theme_classic()

## ----logistic-step-fun---------------------------------------------------
step.fun <- Csnippet("
  double dW = rnorm(0,sqrt(dt));
  N += r*N*(1-N/K)*dt+sigma*N*dW;
")

## ----logistic-pomp1------------------------------------------------------
parus <- pomp(data=parus.dat,time="year",t0=1959,
              rprocess=euler.sim(step.fun=step.fun,delta.t=1/365),
              statenames="N",paramnames=c("r","K","sigma"))

## ----logistic-simul1-----------------------------------------------------
simStates <- simulate(parus,nsim=10,params=c(r=0.2,K=200,sigma=0.5,N.0=200),states=TRUE)

## ----logistic-plot1,echo=FALSE-------------------------------------------
melt(simStates) %>%
  dcast(rep+time~variable) %>%
  ggplot(mapping=aes(x=time,y=N,group=rep,color=factor(rep)))+
  geom_line()+guides(color=FALSE)+
  theme_classic()

## ----logistic-rmeasure---------------------------------------------------
rmeas <- Csnippet("
  P = rpois(N);
")

## ----logistic-pomp2------------------------------------------------------
parus <- pomp(parus,rmeasure=rmeas,statenames="N")

## ----logistic-simul2-----------------------------------------------------
sim <- simulate(parus,params=c(r=0.2,K=200,sigma=0.5,N.0=200),
                nsim=10,obs=TRUE,states=TRUE)

## ----logistic-plot2,echo=FALSE-------------------------------------------
sim %>% melt() %>%
  ggplot(mapping=aes(x=time,y=value,group=rep,color=factor(rep)))+
  geom_line()+
  guides(color=FALSE)+scale_y_sqrt()+
  facet_grid(variable~.,scales="free_y")

sim %>% melt() %>% dcast(rep+time~variable,value.var='value') %>%
  ggplot(mapping=aes(x=N,y=P,color=factor(rep)))+
  geom_point()+scale_x_sqrt()+scale_y_sqrt()+
  coord_equal()+
  guides(color=FALSE)

## ----logistic-dmeasure---------------------------------------------------
dmeas <- Csnippet("
  lik = dpois(P,N,give_log);
")

## ----logistic-pomp3------------------------------------------------------
parus <- pomp(parus,dmeasure=dmeas,statenames="N")

## ----logistic-pfilter----------------------------------------------------
pf <- pfilter(parus,Np=1000,params=c(r=0.2,K=200,sigma=0.5,N.0=200))
logLik(pf)

## ----logistic-skeleton---------------------------------------------------
skel <- Csnippet("
  DN = r*N*(1-N/K);
")

parus <- pomp(parus,skeleton=vectorfield(skel),statenames="N",paramnames=c("r","K"))

## ----logistic-traj1------------------------------------------------------
pars <- parmat(c(r=1,K=200,sigma=0.5,N.0=20),5)
pars["N.0",] <- seq(20,300,length=5)
traj <- trajectory(parus,params=pars,times=seq(1959,1970,by=0.01))

## ----logistic-plot3,echo=FALSE-------------------------------------------
trajectory(parus,params=pars,times=seq(1959,1970,by=0.01),as.data.frame=TRUE) %>%
  ggplot(mapping=aes(x=time,y=N,group=traj,color=traj))+
  guides(color=FALSE)+
  geom_line()

## ----bh-stepfun----------------------------------------------------------
bh.step <- Csnippet("
  double eps = rlnorm(-sigma*sigma/2,sigma);
  N = a*N/(1+b*N)*eps;
")

## ----bh-skeleton---------------------------------------------------------
bh.skel <- Csnippet("
  DN = a*N/(1+b*N);
")

## ----bh-pomp1------------------------------------------------------------
parus.bh <- pomp(parus,rprocess=discrete.time.sim(bh.step,delta.t=1),
                 skeleton=bh.skel,
                 statenames="N",paramnames=c("a","b","sigma"))
parus.bh <- pomp(parus,rprocess=discrete.time.sim(bh.step,delta.t=1),
                 skeleton=map(bh.skel,delta.t=1),
                 statenames="N",paramnames=c("a","b","sigma"))

## ----bh-test-------------------------------------------------------------
coef(parus.bh) <- c(a=1.1,b=5e-4,sigma=0.5,N.0=30)
sim <- simulate(parus.bh)
traj <- trajectory(parus.bh)
pf <- pfilter(parus.bh,Np=1000)

## ----logistic-partrans---------------------------------------------------
logtrans <- Csnippet("
  Tr = log(r);
  TK = log(K);
  Tsigma = log(sigma);
")

exptrans <- Csnippet("
  Tr = exp(r);
  TK = exp(K);
  Tsigma = exp(sigma);
")

parus <- pomp(parus,toEstimationScale=logtrans,
              fromEstimationScale=exptrans,
              paramnames=c("r","K","sigma"))

## ----logistic-partrans-test,include=FALSE--------------------------------
p <- c(r=1,K=200,N.0=200,sigma=0.5)
coef(parus,transform=TRUE) <- partrans(parus,p,dir="inv")
stopifnot(all.equal(p,coef(parus)))

## ----parus-traj-match----------------------------------------------------
tm <- traj.match(parus,start=c(r=1,K=200,N.0=200,sigma=0.5),
                 est=c("r","K"),transform=TRUE)
signif(coef(tm),3)
logLik(tm)

## ----parus-tm-sim1-------------------------------------------------------
coef(tm,"sigma") <- 0
simulate(tm,nsim=10,as.data.frame=TRUE,include.data=TRUE) %>%
  ggplot(aes(x=time,y=P,group=sim,alpha=(sim=="data")))+
  scale_alpha_manual(name="",values=c(`TRUE`=1,`FALSE`=0.2),
                     labels=c(`FALSE`="simulation",`TRUE`="data"))+
  geom_line()

## ----parus-dprior--------------------------------------------------------
parus %<>%
  pomp(dprior=Csnippet("
    lik = dunif(r,0,5,1)+dunif(K,100,600,1)+dunif(sigma,0,2,1);
    lik = (give_log) ? lik : exp(lik);
  "),paramnames=c("r","K","sigma"))

dev.off()

capture.output(show(parus)) -> ignore
