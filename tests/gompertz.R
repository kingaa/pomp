options(digits=3)

library(pomp)
library(magrittr)

pompExample(gompertz)

set.seed(1438408329L)

stopifnot(all.equal(coef(gompertz),partrans(gompertz,coef(gompertz,transform=TRUE),dir="from")))
po <- simulate(gompertz)
stopifnot(round(sum(obs(po)),3)==114.451)
pf <- pfilter(gompertz,Np=1000)
stopifnot(all.equal(round(mean(eff.sample.size(pf)),3),593.45))
tj <- trajectory(gompertz)
stopifnot(all.equal(round(sum(tj["X",,]),3),101))

gompertz %>% window(end=10) %>% simulate(seed=1176423047) -> po

set.seed(376098756)

po %>%
  mif2(Nmif=100,Np=1000,cooling.fraction.50=0.4,cooling.type="geometric",
    rw.sd=rw.sd(sigma=0.02,r=0.02,X.0=ivp(0.05),tau=0.02),
    transform=TRUE) %>%
  continue(Nmif=100) -> mf
replicate(n=10,mf %>% pfilter(Np=3000)) -> pfs
pfs %>% sapply(logLik) %>% logmeanexp(se=TRUE) -> pf.ll.mle

replicate(n=10,po %>% pfilter(Np=3000)) %>%
  sapply(logLik) %>%
  logmeanexp(se=TRUE) -> pf.ll.truth

po %>% init.state(params=coef(po)) %>% as.numeric() -> x0
po %>% coef() %>% as.list() %$% {c(x=log(x0/K))} -> X0
po %>% obs() %>% log() -> y
po %>% coef() %>% as.list() %$% {matrix(c(exp(-r)),1,1)} -> A
po %>% coef() %>% as.list() %$% {matrix(c(sigma*sigma),1,1)} -> Q
po %>% coef() %>% as.list() %$% {matrix(1,1,1)} -> C
po %>% coef() %>% as.list() %$% {tau*tau} -> R
po %>% time() -> t
pomp:::kalmanFilter(t,y,X0,A,Q,C,R) -> kf.truth

mf %>% init.state(params=coef(mf)) %>% as.numeric() -> x0
mf %>% coef() %>% as.list() %$% {c(x=log(x0/K))} -> X0
mf %>% obs() %>% log() -> y
mf %>% coef() %>% as.list() %$% {matrix(c(exp(-r)),1,1)} -> A
mf %>% coef() %>% as.list() %$% {matrix(c(sigma*sigma),1,1)} -> Q
mf %>% coef() %>% as.list() %$% {matrix(1,1,1)} -> C
mf %>% coef() %>% as.list() %$% {tau*tau} -> R
mf %>% time() -> t
pomp:::kalmanFilter(t,y,X0,A,Q,C,R) -> kf.mle

cat("likelihood at truth:",kf.truth$loglik-sum(y),"\n")
cat("pfilter likelihood at truth:",pf.ll.truth,"\n")
cat("likelihood at IF2 mle:",kf.mle$loglik-sum(y),"\n")
cat("pfilter likelihood at IF2 mle:",pf.ll.mle,"\n")

