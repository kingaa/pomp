library(pomp)

pompExample(ou2)

set.seed(64857673L)

obs(window(ou2,end=20,start=15))
obs(window(ou2,end=5),"y1")

fit1.pfilter <- pfilter(ou2,Np=1000)
cat("coefficients at `truth'\n")
print(coef(ou2,c('x1.0','x2.0','alpha.2','alpha.3')),digits=4)
cat("particle filter log likelihood at truth\n")
print(fit1.pfilter$loglik,digits=4)

p.truth <- coef(ou2)
guess2 <- guess1 <- p.truth
guess1[c('x1.0','x2.0','alpha.2','alpha.3')] <- 0.5*guess1[c('x1.0','x2.0','alpha.2','alpha.3')]
guess2[c('x1.0','x2.0','alpha.2','alpha.3')] <- 1.5*guess1[c('x1.0','x2.0','alpha.2','alpha.3')]

mif1 <- mif(ou2,Nmif=30,start=guess1,
            pars=c('alpha.2','alpha.3'),ivps=c('x1.0','x2.0'),
            rw.sd=c(
              x1.0=5,x2.0=5,
              alpha.2=0.1,alpha.3=0.1
              ),
            Np=1000,
            var.factor=1,
            ic.lag=10,
            cooling.type="geometric",
            cooling.fraction=0.95^50,
            max.fail=100
            )

mif2 <- mif(ou2,Nmif=30,start=guess2,
            pars=c('alpha.2','alpha.3'),ivps=c('x1.0','x2.0'),
            transform=TRUE,
            rw.sd=c(
              x1.0=5,x2.0=5,
              alpha.2=0.1,alpha.3=0.1
              ),
            Np=1000,
            var.factor=1,
            ic.lag=10,
            cooling.type="geometric",
            cooling.fraction=0.95^50,
            max.fail=100
            )

pdf(file="ou2-mif.pdf")
plot(mif1)
plot(mif2)
try(plot(mif1,mif2))
plot(c(mif1,mif2))
coef(c(mif1,mif2))
dev.off()

set.seed(33848585L)

try(
    mif(
        ou2,
        Nmif=1,
        pars=c("alpha.1","alpha.4","x1.0"),
        ivps=c("x1.0","x2.0"),
        rw.sd=c(x1.0=5,x2.0=5,alpha.1=0.1,alpha.4=0.2,alpha.3=0),
        Np=100,cooling.type="geometric",cooling.fraction=0.95^50,
        ic.lag=10,var.factor=1
        )
    )

try(
    mif(
        ou2,
        Nmif=1,
        pars=c("alpha.1","alpha.4"),
        ivps=c("x1.0","x2.0"),
        rw.sd=c(x1.0=5,x2.0=5,alpha.1=0,alpha.4=0.2,alpha.3=0),
        Np=100,
        cooling.type="geometric",cooling.fraction=0.95^50,
        ic.lag=10,var.factor=1
        )
    )

try(
    mif(
        ou2,
        Nmif=1,
        ivps=c("x1.0","x2.0"),
        rw.sd=c(x1.0=5,x2.0=5,alpha.1=0.1,alpha.4=0.2,alpha.3=0),
        Np=100,ic.lag=10,var.factor=1
        )
    )

try(
    mif(
        ou2,
        Nmif=1,
        ivps=c("x1.0","x2.0"),
        rw.sd=c(x1.0=5,x2.0=5,alpha.1=0,alpha.4=0.2,alpha.3=0),
        Np=-10,cooling.type="geometric",cooling.fraction=0.95^50,
        ic.lag=10,var.factor=1
        )
    )

try(
    mif(
        ou2,
        Nmif=-3,
        ivps=c("x1.0","x2.0"),
        rw.sd=c(x1.0=5,x2.0=5,alpha.1=0,alpha.4=0.2,alpha.3=0),
        Np=11.6,cooling.type="geometric",cooling.fraction=0.95^50,
        ic.lag=10,var.factor=1
        )
    )

try(
    mif(
        ou2,
        Nmif=2,
        start=c(alpha.1=0.9,alpha.2=0,alpha.3=0,alpha.4=-Inf,sigma.1=1,sigma.2=0,sigma.3=2,tau=1,x1.0=50,x2.0=-50),
        ivps=c("x1.0","x2.0"),
        rw.sd=c(x1.0=5,x2.0=5,alpha.1=0,alpha.4=0.2,alpha.3=0),
        Np=11,cooling.type="geometric",cooling.fraction=0.95^50,
        ic.lag=10,var.factor=1
        )
    )

try(
    mif(
        ou2,
        Nmif=2,
        start=c(alpha.1=0.9,alpha.2=0,alpha.3=0,alpha.4=0.99,sigma.1=1,sigma.2=0,sigma.3=2,tau=1,x1.0=50,x2.0=NaN),
        ivps=c("x1.0","x2.0"),
        rw.sd=c(x1.0=5,x2.0=5,alpha.1=0,alpha.4=0.2,alpha.3=0),
        Np=11,cooling.type="geometric",cooling.fraction=0.95^50,
        ic.lag=10,var.factor=1
        )
    )

fit <- mif(
           ou2,
           Nmif=0,
           pars=c("alpha.2","alpha.3"),
           ivps=c("x1.0","x2.0"),
           rw.sd=c(x1.0=5,x2.0=5,alpha.2=0.1,alpha.3=0.2,alpha.3=0),
           Np=100,cooling.type="geometric",cooling.fraction=0.95^50,
           ic.lag=10,var.factor=1
           )
fit <- mif(
           fit,
           Nmif=2,
           ivps=c("x1.0","x2.0"),
           rw.sd=c(x1.0=5,x2.0=5,alpha.2=0.1,alpha.3=0.2),
           Np=function(k)if(k<10) 2000 else 500,
           cooling.type="geometric",cooling.fraction=0.95^50,
           ic.lag=10,var.factor=1
           )
fit <- continue(fit)
fit <- continue(fit,Nmif=2)
ff <- pfilter(fit,pred.mean=T,filter.mean=T,pred.var=T,max.fail=100,verbose=T)
ff <- pfilter(ff)
fit <- mif(fit,rw.sd=c(x1.0=5,x2.0=5,alpha.2=0.1,alpha.3=0.1))
fit <- continue(fit,Nmif=2,ivps=c("x1.0"),pars=c("alpha.2"))
s <- coef(fit)
s[2] <- 0.01
fit <- mif(fit,Nmif=3,start=s)
fit <- mif(ou2,Nmif=3,rw.sd=c(alpha.2=0.1,alpha.3=0.1),Np=1000,cooling.type="geometric",cooling.fraction=0.98^50,var.factor=1,ic.lag=2)
fit <- continue(fit,Nmif=2,Np=2000)
fit <- continue(fit,ivps=c("x1.0"),rw.sd=c(alpha.2=0.1,alpha.3=0.1,x1.0=5,x2.0=5),Nmif=2)
ff <- pfilter(fit)
fit <- mif(
           ff,
           Nmif=2,
           ivps=c("x1.0","x2.0"),
           rw.sd=c(x1.0=5,x2.0=5,alpha.2=0.1,alpha.3=0.2),
           cooling.fraction=0.95^50,ic.lag=10,var.factor=1
           )

pp <- particles(fit,Np=10,center=coef(fit),sd=abs(0.1*coef(fit)))
fit <- mif(
           fit,
           Nmif=10,
           Np=1000,
           particles=function(Np,center,sd,...){
             matrix(
                    data=runif(
                      n=Np*length(center),
                      min=center-sd,
                      max=center+sd
                      ),
                    nrow=length(center),
                    ncol=Np,
                    dimnames=list(names(center),NULL)
                    )
           }
           )
pp <- particles(fit,Np=10,center=coef(fit),sd=abs(0.1*coef(fit)))

