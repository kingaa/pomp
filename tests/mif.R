library(pomp)

pompExample(ou2)

png(file="mif-%02d.png",res=100)

set.seed(64857673L)
options(digits=3)

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
            ivps=c('x1.0','x2.0'),
            rw.sd=c(
              x1.0=5,x2.0=5,
              alpha.2=0.1,alpha.3=0.1
              ),
            Np=1000,
            var.factor=1,
            ic.lag=10,
            cooling.type="hyperbolic",
            cooling.fraction=0.95^50,
            max.fail=100
            )

capture.output(
    mif2 <- mif(ou2,Nmif=30,start=guess2,
                ivps=c('x1.0','x2.0'),
                transform=TRUE,
                verbose=TRUE,
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
    ) -> out
stopifnot(length(out)==630)
stopifnot(sum(grepl("pfilter timestep",out))==600)
stopifnot(sum(grepl("mif iteration",out))==30)

plot(mif1)
plot(mif12 <- c(mif1,mif2))
mif12 <- c(mif12)
mif1a <- c(mif1)
coef(mif2)
dim(coef(mif12))
dim(coef(c(mif12,mif2)))
dim(coef(c(mif1,mif12)))
dim(coef(c(mif1a,mif12)))
dim(coef(mif12[2]))
dim(conv.rec(mif12))

set.seed(33848585L)

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
try(continue(fit,Np=function(k)if(k<10) c(20,30) else 500))

fit <- continue(fit)
fit <- continue(fit,Nmif=2)
capture.output(
    ff <- pfilter(fit,pred.mean=T,filter.mean=T,pred.var=T,max.fail=100,verbose=TRUE)) -> out
stopifnot(length(out)==20)
stopifnot(sum(grepl("pfilter timestep",out))==20)
ff <- pfilter(ff)
fit <- mif(fit,rw.sd=c(x1.0=5,x2.0=5,alpha.2=0.1,alpha.3=0.1))
fit <- continue(fit,Nmif=2,ivps=c("x1.0"))
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
try(mif(fit,Nmif=2,method="mif2"))
fit <- continue(fit,Nmif=1,ic.lag=10000)
fit <- continue(fit,Nmif=1,ic.lag=2,rw.sd=c(x1.0=5,x2.0=5))
fit <- continue(fit,Nmif=3,ic.lag=100,method="unweighted")
try(mif(ff,Nmif=2,ivps=c("x1.0","x2.0"),
        rw.sd=c(x1.0=5,x2.0=5,alpha.2=0.1,alpha.3=0.2),
        cooling.fraction=0.95^50,var.factor=1))
coef(fit,"alpha.4") <- 0
try(fit <- continue(fit,Np=2,Nmif=3,max.fail=1))
coef(fit,"tau") <- NaN
try(fit <- continue(fit,Np=100,Nmif=2))
coef(fit,"alpha.2") <- -Inf
try(fit <- continue(fit,Np=100,Nmif=2))

pompExample(gompertz)

coef(gompertz,"K") <- -1
try(mif(gompertz,Np=1000,rw.sd=c(K=0.1,r=0.1),cooling.fraction.50=0.5))
try(mif(gompertz,Np=1000,rw.sd=rw.sd(K=0.1,r=0.1),cooling.fraction.50=0.5))

pomp(gompertz,
     toEstimationScale=function (params,...){
         params["r"] <- log(params["r"])
         params
     },
     fromEstimationScale=function (params,...){
         params["r"] <- exp(params["r"])
         params
     }) -> po

try(mif(po,Nmif=3,Np=1000,rw.sd=c(K=5,r=5),cooling.fraction.50=0.5))
try(mif(po,Nmif=3,Np=1000,rw.sd=c(K=5,r=Inf),cooling.fraction.50=0.5))

coef(gompertz,"K") <- 10
try(gb <- mif(gompertz,Nmif=1,Np=1,rw.sd=c(K=0.1,r=0.1),
              cooling.fraction.50=0.2,transform=TRUE))

dev.off()
