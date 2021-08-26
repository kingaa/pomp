options(digits=3)

library(pomp)
library(magrittr)
library(dplyr)

set.seed(376098756)

gompertz() %>% window(end=10) %>% simulate(seed=1176423047) -> po

po %>%
  mif2(Nmif=100,Np=1000,cooling.fraction.50=0.4,cooling.type="geometric",
    rw.sd=rw.sd(sigma=0.02,r=0.02,X_0=ivp(0.05),tau=0.02)) %>%
  continue(Nmif=100) -> mf
replicate(n=10,mf %>% pfilter(Np=3000)) -> pfs
pfs %>% sapply(logLik) %>% logmeanexp(se=TRUE) -> pf.ll.mle

replicate(n=10,po %>% pfilter(Np=3000)) %>%
  sapply(logLik) %>%
  logmeanexp(se=TRUE) -> pf.ll.truth

po %>%
  as.data.frame() %>%
  mutate(logY=log(Y)) %>%
  select(time,logY) %>%
  pomp(t0=0,times="time") -> logpo
  
po %>% rinit(params=coef(po)) %>% as.numeric() -> x0
po %>% coef() %>% as.list() %$% {c(x=log(x0/K))} -> X0
po %>% obs() %>% log() -> y
po %>% coef() %>% as.list() %$% {matrix(c(exp(-r)),1,1)} -> A
po %>% coef() %>% as.list() %$% {matrix(c(sigma*sigma),1,1)} -> Q
po %>% coef() %>% as.list() %$% {matrix(1,1,1)} -> C
po %>% coef() %>% as.list() %$% {tau*tau} -> R
kalmanFilter(logpo,X0=X0,A=A,Q=Q,C=C,R=R) -> kf.truth

mf %>% rinit(params=coef(mf)) %>% as.numeric() -> x0
mf %>% coef() %>% as.list() %$% {c(x=log(x0/K))} -> X0
mf %>% obs() %>% log() -> y
mf %>% coef() %>% as.list() %$% {matrix(c(exp(-r)),1,1)} -> A
mf %>% coef() %>% as.list() %$% {matrix(c(sigma*sigma),1,1)} -> Q
mf %>% coef() %>% as.list() %$% {matrix(1,1,1)} -> C
mf %>% coef() %>% as.list() %$% {tau*tau} -> R
kalmanFilter(logpo,X0=X0,A=A,Q=Q,C=C,R=R) -> kf.mle

cat("likelihood at truth:",kf.truth$logLik-sum(y),"\n")
cat("pfilter likelihood at truth:",pf.ll.truth,"\n")
cat("likelihood at IF2 mle:",kf.mle$logLik-sum(y),"\n")
cat("pfilter likelihood at IF2 mle:",pf.ll.mle,"\n")
