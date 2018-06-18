library(pomp)

pompExample(ou2)

png(file="mif2-%02d.png",res=100)

set.seed(64857673L)
options(digits=3)

guess2 <- guess1 <- coef(ou2)
guess1[c('x1.0','x2.0','alpha.2','alpha.3')] <- 0.5*guess1[c('x1.0','x2.0','alpha.2','alpha.3')]
guess2[c('x1.0','x2.0','alpha.2','alpha.3')] <- 1.5*guess1[c('x1.0','x2.0','alpha.2','alpha.3')]

try(mif2(ou2,Nmif=10,start=guess2,
         cooling.type="hyperbolic",cooling.fraction.50=0.05,
         rw.sd=rw.sd(
             x1.0=ivp(0.5),x2.0=ivp(0.5),
             alpha.2=0.1,alpha.3=0.1)) -> ignore)
try(mif2(ou2,Nmif=1,start=guess2,Np=rep(10,102),
         cooling.type="hyperbolic",cooling.fraction.50=0.05,
         rw.sd=rw.sd(
             x1.0=ivp(0.5),x2.0=ivp(0.5),
             alpha.2=0.1,alpha.3=0.1)) -> ignore)
try(mif2(ou2,Nmif=1,start=guess2,Np=seq_len(101),
         cooling.type="hyperbolic",cooling.fraction.50=0.05,
         rw.sd=rw.sd(
             x1.0=ivp(0.5),x2.0=ivp(0.5),
             alpha.2=0.1,alpha.3=0.1)) -> ignore)
try(mif2(ou2,Nmif=0,start=guess2,Np=1000,
         cooling.type="hyperbolic",cooling.fraction.50=0.05,
         rw.sd=rw.sd(
             x1.0=ivp(0.5),x2.0=ivp(0.5),
             alpha.2=0.1,alpha.3=0.1)) -> ignore)
try(mif2(ou2,Nmif=1,Np=1000,rw.sd=rw.sd(ivp(0.1))))
try(mif2(ou2,Nmif=1,Np=10,rw.sd=rw.sd(alpha.2=10),
  cooling.fraction.50=1,max.fail=0))
try(mif2(ou2,Nmif=1,Np=100,start=numeric(0)))
try(mif2(ou2,Nmif=1,Np=100,start=unname(coef(ou2))))
try(mif2(ou2,Nmif=1,Np="100"))
try(mif2(ou2,Nmif=1,Np=10,cooling.fraction.50=1))
try(mif2(ou2,Nmif=1,Np=10,rw.sd=rw.sd(alpha.2=10),
  cooling.fraction.50=2,max.fail=0))

m1 <- mif2(ou2,Nmif=50,start=guess1,Np=1000,
           cooling.type="hyperbolic",cooling.fraction.50=0.05,
           rw.sd=rw.sd(
               x1.0=ivp(0.5),x2.0=ivp(0.5),
               alpha.2=0.1,alpha.3=0.1))
m2 <- mif2(ou2,Nmif=30,start=guess2,Np=1000,
           cooling.type="hyperbolic",cooling.fraction.50=0.05,
           rw.sd=rw.sd(
               x1.0=ivp(0.5),x2.0=ivp(0.5),
               alpha.2=0.1,alpha.3=0.1))
m2 <- continue(m2,Nmif=20)

try(conv.rec(m2,c("bob","nancy")))

freeze(continue(m2,Nmif=2,.paramMatrix=m2@paramMatrix),seed=595996) -> m2a
freeze(continue(m2,Nmif=2),seed=595996) -> m2b
stopifnot(all.equal(coef(m2a),coef(m2b)))

try(continue(m2,Nmif=10,rw.sd=rw.sd(beta=0.1,alpha.2=0.1)))
try(continue(m2,Nmif=10,rw.sd=rw.sd(alpha.2=rep(0.1,3))))

plot(m1,y=NA)
dim(coef(c(m1)))
plot(m12 <- c(m1,m2),y=33)
coef(c(m12))
dim(coef(m12))
dim(coef(c(m12,m12)))
dim(coef(c(m12,m1)))
dim(coef(c(m1,m12)))
dim(coef(m12[2]))
sapply(conv.rec(m12),dim)
try(c(m1,ou2))
try(c(c(m1,m2),ou2))

rbind(mle1=c(coef(m1),loglik=logLik(pfilter(m1,Np=1000))),
      mle2=c(coef(m2),loglik=logLik(pfilter(m1,Np=1000))),
      truth=c(coef(ou2),loglik=logLik(pfilter(m1,Np=1000))))

m3 <- mif2(ou2,Nmif=3,start=guess1,Np=200,
           cooling.fraction=0.2,
           rw.sd=rw.sd(
               x1.0=c(0.5,rep(0.2,99)),
               x2.0=ivp(0.5),
               alpha.2=ifelse(time==1,0.2,0.1),
               alpha.3=0.2*(time<10)))

m4 <- mif2(ou2,Nmif=3,start=guess1,
           Np=function(k)if(k<20) 200 else 100,
           cooling.fraction=0.2,
           rw.sd=rw.sd(
               x1.0=c(0.5,rep(0.2,99)),
               x2.0=ivp(0.5),
               alpha.2=ifelse(time==1,0.2,0.1),
               alpha.3=0.2*(time<10)))

m4 <- mif2(m4,Nmif=10,Np=c(rep(200,20),rep(100,80),200))
m4 <- continue(m4,Nmif=10,cooling.fraction=0.1)
try(continue(m4,Np=function(k)if(k<10) "B" else 500))
try(continue(m4,Np=function(k)if(k<10) -30 else 500))

pf <- pfilter(m4)
half <- 0.5
capture.output(
    m4 <- mif2(pf,Nmif=10,
               cooling.fraction=0.2,
               verbose=TRUE,
               rw.sd=rw.sd(
                   x1.0=c(half,rep(0.2,99)),
                   x2.0=ivp(half),
                   alpha.2=ifelse(time==1,0.2,0.1),
                   alpha.3=0.2*(time<10)))
) -> out
stopifnot(length(out)==210)
stopifnot(sum(grepl("mif2 iteration",out))==10)
stopifnot(sum(grepl("mif2 pfilter",out))==200)

m5 <- m4
coef(m5,"alpha.2") <- -Inf
try(mif2(m5,Np=100,Nmif=2))

library(ggplot2)
library(reshape2)
library(magrittr)

m4 %>% conv.rec() %>% melt() %>%
    ggplot(aes(x=iteration,y=value,color=variable))+
    geom_line()+
    facet_wrap(~variable,scales='free_y',ncol=2)

m4 %>% conv.rec(c("alpha.2","alpha.3","loglik")) %>% melt() %>%
    ggplot(aes(x=iteration,y=value,color=variable))+
    geom_line()+
    facet_wrap(~variable,scales="free_y",ncol=1)

pompExample(gompertz)

coef(gompertz,"K") <- -1
try(mif2(gompertz,Np=1000,rw.sd=rw.sd(K=0.1,r=0.1),cooling.fraction.50=0.5))

pomp(gompertz,
     toEstimationScale=function (params,...){
         params["r"] <- log(params["r"])
         params
     },
     fromEstimationScale=function (params,...){
         params["r"] <- exp(params["r"])
         params
     }) -> po

try(mif2(po,Nmif=3,Np=1000,rw.sd=rw.sd(K=5,r=Inf),cooling.fraction.50=0.5))

coef(gompertz,"K") <- 10
try(gb <- mif2(gompertz,Nmif=1,Np=1,rw.sd=c(K=0.1,r=0.1),
               cooling.fraction.50=0.2,transform=TRUE))

mif2(m4,Nmif=2,cooling.type="hyperbolic",cooling.fraction.50=1) -> m6

ep <- 0.1
m5 <- mif2(ou2,Nmif=2,start=guess2,Np=1000,
           cooling.type="hyperbolic",cooling.fraction.50=0.05,
           rw.sd=rw.sd(x1.0=ivp(ep),x2.0=ivp(ep),alpha.2=ep,alpha.3=ep))
stopifnot(is(m5,"mif2d.pomp"))

f <- function () {
    ep <- 0.2
    m <- mif2(ou2,Nmif=2,start=guess2,Np=1000,
              cooling.type="hyperbolic",cooling.fraction.50=0.05,
              rw.sd=rw.sd(x1.0=ivp(ep),alpha.2=ep))
    as.numeric(m@rw.sd[,1:2])
}

stopifnot(all.equal(f(),c(0.2,0.2,0.0,0.2)))

rr <- rw.sd(x1.0=ivp(ep),alpha.2=ep)
f <- function () {
    m <- mif2(ou2,Nmif=2,start=guess2,Np=1000,
              cooling.type="hyperbolic",cooling.fraction.50=0.05,
              rw.sd=rr)
    as.numeric(m@rw.sd[,1:2])
}

stopifnot(all.equal(f(),c(0.1,0.1,0.0,0.1)))

m7 <- m1
coef(m7,"tau") <- 0
m7 <- mif2(m7,Nmif=2)
m7 <- mif2(m7)

dev.off()
