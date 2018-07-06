library(pomp)

pompExample(ou2)

set.seed(9994847L)

png(filename="pfilter-%02d.png",res=100)

pf <- pfilter(ou2,Np=1000)
print(coef(ou2,c('x1.0','x2.0','alpha.1','alpha.4')),digits=4)
cat("particle filter log likelihood at truth\n")
print(pf$loglik,digits=4)

try(pfilter(ou2,Np=1000,params=c()))

try(pfilter(ou2,Np=c(1000,10)))
try(pfilter(ou2,Np=-100))
try(pfilter(ou2,Np="bob"))
try(pfilter(ou2,Np=100,params=matrix(0,nrow=3,ncol=100)))
try(pfilter(ou2,Np=1,max.fail=0))
try(replicate(n=10,pfilter(ou2,Np=function(k)if(k<10) c(100,200) else 500)))
pf <- replicate(n=10,pfilter(ou2,Np=function(k)if(k<10) 10000 else 500))
pf.ll <- sapply(pf,logLik)
ll.est <- log(mean(exp(pf.ll-mean(pf.ll))))+mean(pf.ll)
ll.se <- sd(exp(pf.ll-mean(pf.ll)))/exp(ll.est-mean(pf.ll))/sqrt(length(pf))
print(round(c(loglik=ll.est,loglik.se=ll.se),digits=2))

pompExample(euler.sir)
pf <- pfilter(euler.sir,Np=100)
print(coef(pf))
print(pf$loglik,digits=4)

coef(pf,"rho") <- -1
try(pfilter(pf))

p <- coef(euler.sir)
euler.sir@params <- numeric(0)
p["iota"] <- 1
pf <- pfilter(euler.sir,params=p,Np=100,filter.mean=TRUE)
print(coef(pf))
print(logLik(pf),digits=4)
plot(cond.loglik~time,data=as(pf,"data.frame"),type='l')
plot(ess~time,data=as.data.frame(pf),type='l')
plot(filter.mean.I~time,data=as(pf,"data.frame"),type='l')

pompExample(gompertz)
pfilter(gompertz,params=parmat(coef(gompertz),100)) -> pf
pfilter(gompertz,params=parmat(coef(gompertz),100),Np=100) -> pf
names(as.data.frame(pfilter(gompertz,Np=100,pred.mean=TRUE,pred.var=TRUE)))
pfilter(gompertz,params=as.list(coef(gompertz)),Np=100) -> pf
try(pfilter(gompertz) -> pf)
try(pfilter(gompertz,params=parmat(coef(gompertz),100),Np=1000) -> pf)
coef(gompertz,"sigma") <- Inf
try(pfilter(gompertz,Np=1000,pred.var=TRUE))

dev.off()
