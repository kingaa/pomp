## ----prelims,echo=F,cache=F----------------------------------------------
require(pomp)
options(
  keep.source=TRUE,
  stringsAsFactors=FALSE,
  encoding="UTF-8",
  scipen=5
  )
set.seed(588665L)

## ----eval=FALSE----------------------------------------------------------
## install.packages("pomp",repos="http://R-Forge.R-Project.org")

## ----eval=FALSE----------------------------------------------------------
## cat("#include <R.h>
##     void hello (void) {
##     Rprintf(\"hello world!\\n\");
##     }",file="hello.c")
## system("R CMD SHLIB hello.c")
## dyn.load(paste0("hello",.Platform$dynlib.ext))
## .C("hello",PACKAGE="hello")

## ----load-parus-data-----------------------------------------------------
parus.dat <- read.csv(text="
                      year,pop
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

## ----logistic-step-fun---------------------------------------------------
step.fun <- Csnippet("
  double dW = rnorm(0,sqrt(dt));
  N += r*N*(1-N/K)*dt+sigma*N*dW;
")

## ----logistic-pomp1,cache=FALSE------------------------------------------
parus <- pomp(data=parus.dat,time="year",t0=1959,
              rprocess=euler.sim(step.fun=step.fun,delta.t=1/365),
              statenames="N",paramnames=c("r","K","sigma"))

## ----logistic-simul1-----------------------------------------------------
simStates <- simulate(parus,nsim=10,params=c(r=0.2,K=200,sigma=0.5,N.0=200),states=TRUE)

## ----logistic-rmeasure---------------------------------------------------
rmeas <- Csnippet("
  pop = rpois(phi*N);
")

## ----logistic-pomp2,cache=FALSE------------------------------------------
parus <- pomp(parus,rmeasure=rmeas,paramnames="phi",statenames="N")

## ----logistic-simul2-----------------------------------------------------
sim <- simulate(parus,params=c(r=0.2,K=200,phi=0.5,sigma=0.5,N.0=200),
                nsim=10,obs=TRUE,states=TRUE)

## ----logistic-dmeasure---------------------------------------------------
dmeas <- Csnippet("
  lik = dpois(pop,phi*N,give_log);
")

## ----logistic-pomp3,cache=FALSE------------------------------------------
parus <- pomp(parus,dmeasure=dmeas,paramnames="phi",statenames="N")

## ----logistic-pfilter----------------------------------------------------
pf <- pfilter(parus,Np=1000,params=c(r=0.2,K=200,phi=0.5,sigma=0.5,N.0=200))
logLik(pf)

## ----logistic-skeleton,cache=FALSE---------------------------------------
skel <- Csnippet("
  DN = r*N*(1-N/K);
")

parus <- pomp(parus,skeleton=skel,skeleton.type="vectorfield",statenames="N",paramnames=c("r","K"))

## ----logistic-traj1------------------------------------------------------
pars <- parmat(c(r=1,K=200,phi=0.5,sigma=0.5,N.0=20),5)
pars["N.0",] <- seq(20,300,length=5)
traj <- trajectory(parus,params=pars,times=seq(1959,1970,by=0.01))

## ----bh-stepfun----------------------------------------------------------
bh.step <- Csnippet("
  double eps = rlnorm(-sigma*sigma/2,sigma);
  N = a*N/(1+b*N)*eps;
")

## ----bh-skeleton---------------------------------------------------------
bh.skel <- Csnippet("
  DN = a*N/(1+b*N);
")

## ----bh-pomp1,cache=FALSE------------------------------------------------
parus.bh <- pomp(parus,rprocess=discrete.time.sim(bh.step,delta.t=1),skeleton=bh.skel,skeleton.type="map",skelmap.delta.t=1,statenames="N",paramnames=c("a","b","sigma"))

## ----bh-test-------------------------------------------------------------
coef(parus.bh) <- c(a=1.1,b=5e-4,sigma=0.5,N.0=30,phi=1)
sim <- simulate(parus.bh)
traj <- trajectory(parus.bh)
pf <- pfilter(parus.bh,Np=1000)

## ----logistic-partrans,cache=FALSE---------------------------------------
logtrans <- Csnippet("
  Tr = log(r);
  TK = log(K);
")

exptrans <- Csnippet("
  Tr = exp(r);
  TK = exp(K);
")

parus <- pomp(parus,toEstimationScale=logtrans,
              fromEstimationScale=exptrans,
              paramnames=c("r","K"))

## ----logistic-partrans-test,include=FALSE,cache=FALSE--------------------
p <- c(r=1,K=200,phi=1,N.0=200,sigma=0.5)
coef(parus,transform=TRUE) <- partrans(parus,p,dir="inv")
stopifnot(all.equal(p,coef(parus)))

## ----parus-traj-match,cache=FALSE----------------------------------------
tm <- traj.match(parus,start=c(r=1,K=200,phi=1,N.0=200,sigma=0.5),
                 est=c("r","K","phi"),transform=TRUE)
signif(coef(tm),3)
logLik(tm)

## ----parus-tm-sim1,cache=FALSE-------------------------------------------
coef(tm,"sigma") <- 0
sim1 <- simulate(tm,nsim=10,as.data.frame=TRUE,include.data=TRUE)

## ----parus-mif,cache=TRUE------------------------------------------------
mf <- mif(parus,Nmif=50,Np=1000,method="mif2",cooling.fraction=0.8,
          rw.sd=c(r=0.02,K=0.02,phi=0.02,sigma=0.02),transform=TRUE)
mf <- mif(mf)
mle <- coef(mf); mle
logmeanexp(replicate(5,logLik(pfilter(mf))),se=TRUE)
sim2 <- simulate(mf,nsim=10,as.data.frame=TRUE,include.data=TRUE)
