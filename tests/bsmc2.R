png(filename="bsmc2-%02d.png",res=100)
options(digits=2)

library(pomp)
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(tidyverse))

pompExample(gompertz)
set.seed(398585L)

time(gompertz) <- 1:10

Np <- 20000

try(bsmc2(gompertz,Np=1000,smooth=0.01,est=c("r","K"),tol=1e-6))

smc <- bsmc2(gompertz,rprior=Csnippet("
              K = runif(0.1,1);
              r = rlnorm(log(0.2),1);
              sigma = rlnorm(log(0.1),0.5);"),
             paramnames=c("r","K","sigma"),Np=1000,smooth=0.01,
             est=c("r","K","sigma"),tol=1e-6,transform=TRUE)

plot(smc,y=NA)
plot(smc,pars=c("r","K"),thin=50)
try(plot(smc,pars="bob"))
plot(smc,pars="K")
try(plot(smc,pars=NULL))
range(eff.sample.size(smc))
logEvidence(smc)
range(cond.logEvidence(smc))
stopifnot(sum(cond.logEvidence(smc))==logEvidence(smc))

full_join(
  smc %>% prior_samples() %>% melt(value.name="prior"),
  smc %>% posterior_samples() %>% melt(value.name="posterior"),
  by=c("variable","rep")
) %>%
  gather(dist,value,prior,posterior) %>%
  ggplot(aes(x=value,group=dist,fill=dist,color=dist))+
  geom_density(alpha=0.5)+
  facet_grid(~variable,scales="free_x")+
  theme_bw()

po <- smc
coef(po) <- NULL
try(bsmc2(smc))
try(bsmc2(po))
try(bsmc2(po,params=NULL))
try(bsmc2(po,params="yes"))
try(bsmc2(po,params=list()))
try(bsmc2(po,params=c(1,2,3)))
theta <- coef(smc)
try(bsmc2(po,params=as.list(theta)))
try(bsmc2(po,params=as.list(theta),est="bob"))
try(bsmc2(po,params=as.list(theta),est="r"))
pp <- parmat(theta,100)
capture.output(invisible(bsmc2(po,params=pp,est="r",Np=1,verbose=TRUE))) -> out
stopifnot(sum(grepl("prior.mean",out))==10)
try(bsmc2(po,params=pp,est="r",smooth=2))
try(bsmc2(po,params=pp,est="r",smooth=NA))
try(bsmc2(po,params=pp,est="r",smooth=NULL))
try(bsmc2(po,params=pp,est="r",smooth=list(3)))
try(bsmc2(po,params=pp,est="r",smooth=c(0.1,0.5)))
rownames(pp) <- NULL
try(bsmc2(po,params=pp,est="r"))
try(bsmc2(smc,Np=-10,est="r"))
try(bsmc2(smc,est="r",Np=100,tol=c(3,5)))
try(bsmc2(smc,est="r",Np=100,max.fail=-1,tol=10))
invisible(bsmc2(smc,est="r",Np=100,max.fail=Inf,tol=10))
try(bsmc2(smc,params=theta,est="r",Np=100,
          dmeasure=Csnippet("error(\"whoa nelly!\");")))

theta <- coef(gompertz)
theta["K"] <- 1
try(capture.output(bsmc2(po,Np=2,params=theta,tol=1,max.fail=1,verbose=TRUE)) -> out)


dev.off()
