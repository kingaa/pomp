library(pomp)

pompExample(ou2)

pdf(file="ou2-mif2.pdf")

set.seed(64857673L)

guess2 <- guess1 <- coef(ou2)
guess1[c('x1.0','x2.0','alpha.2','alpha.3')] <- 0.5*guess1[c('x1.0','x2.0','alpha.2','alpha.3')]
guess2[c('x1.0','x2.0','alpha.2','alpha.3')] <- 1.5*guess1[c('x1.0','x2.0','alpha.2','alpha.3')]

m1 <- mif2(ou2,Nmif=100,start=guess1,Np=1000,
           cooling.type="hyperbolic",cooling.fraction.50=0.05,
           rw.sd=rw.sd(
             x1.0=ivp(0.5),x2.0=ivp(0.5),
             alpha.2=0.1,alpha.3=0.1))

m2 <- mif2(ou2,Nmif=50,start=guess2,Np=1000,
           cooling.type="hyperbolic",cooling.fraction.50=0.05,
           rw.sd=rw.sd(
             x1.0=ivp(0.5),x2.0=ivp(0.5),
             alpha.2=0.1,alpha.3=0.1))
m2 <- continue(m2,Nmif=50)

plot(c(m1,m2))
coef(c(m1,m2))

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

m4 <- mif2(m4,Nmif=2,Np=c(rep(200,20),rep(100,80),200))
m4 <- continue(m4,Nmif=2,cooling.fraction=0.1)
pf <- pfilter(m4)
half <- 0.5
m4 <- mif2(pf,Nmif=2,
           cooling.fraction=0.2,
           rw.sd=rw.sd(
             x1.0=c(half,rep(0.2,99)),
             x2.0=ivp(half),
             alpha.2=ifelse(time==1,0.2,0.1),
             alpha.3=0.2*(time<10)))

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

dev.off()
