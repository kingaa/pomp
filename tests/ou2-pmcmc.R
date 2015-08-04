library(pomp)

set.seed(1178744046L)

pompExample(ou2)

dprior.ou2 <- function (params, log, ...) {
  f <- sum(dunif(params,min=coef(ou2)-1,max=coef(ou2)+1,log=TRUE))
  if (log) f else exp(f)
}

pdf(file="ou2-pmcmc.pdf")

f1 <- pmcmc(
            pomp(ou2,dprior=dprior.ou2),
            Nmcmc=20,
            proposal=mvn.diag.rw(c(alpha.2=0.001,alpha.3=0.001)),
            Np=100,
            verbose=FALSE
            )
f1 <- continue(f1,Nmcmc=20)
plot(f1)

ff <- pfilter(f1)
f2 <- pmcmc(
            ff,
            Nmcmc=20,
            rw.sd=c(alpha.2=0.01,alpha.3=0.01),
            verbose=FALSE
            )

f3 <- pmcmc(f2)
f4 <- continue(f3,Nmcmc=20)
filter.traj(f4)[,40,95:100]

plot(c(f2,f3))

try(ff <- c(f3,f4))

if (Sys.getenv("POMP_FULL_TESTS")=="yes") {
  f2a <- pmcmc(f1,Nmcmc=1000,Np=100,verbose=FALSE)
  plot(f2a)
  runs <- rle(as.numeric(conv.rec(f2a,'loglik')))$lengths
  plot(sort(runs))
  acf(conv.rec(f2a,c("alpha.2","alpha.3")))
}

f5 <- pmcmc(
            pomp(ou2,
                 dprior=function (params, log, ...) {
                   f <- sum(dnorm(params,mean=coef(ou2),sd=1,log=TRUE))
                   if (log) f else exp(f)
                 }
                 ),
            Nmcmc=20,
            proposal=mvn.diag.rw(c(alpha.2=0.001,alpha.3=0.001)),
            Np=100,
            verbose=FALSE
            )
f6 <- continue(f5,Nmcmc=20)
plot(f6)

ff <- c(f4,f6)
plot(ff)
plot(conv.rec(ff,c("alpha.2","alpha.3","loglik")))

ff <- c(f2,f3)

try(ff <- c(ff,f4,f6))
try(ff <- c(f4,ou2))
try(ff <- c(ff,ou2))

plot(ff <- c(ff,f5))
plot(conv.rec(c(f2,ff),c("alpha.2","alpha.3")))
plot(conv.rec(ff[2],c("alpha.2")))
plot(conv.rec(ff[2:3],c("alpha.3")))
plot(window(conv.rec(ff[2:3],c("alpha.3")),thin=3,start=2))
plot(conv.rec(ff[[3]],c("alpha.3")))

sig <- array(data=c(0.1,-0.1,0,0.01),
             dim=c(2,2),
             dimnames=list(
               c("alpha.2","alpha.3"),
               c("alpha.2","alpha.3")))
sig <- crossprod(sig)

f7 <- pmcmc(
            pomp(ou2,
                 dprior=function (params, log, ...) {
                   f <- sum(dnorm(params,mean=coef(ou2),sd=1,log=TRUE))
                   if (log) f else exp(f)
                 }
                 ),
            Nmcmc=30,
            proposal=mvn.rw(sig),
            Np=100,
            verbose=FALSE
            )
plot(f7)
filter.traj(f7,"x1")[1,30,1:5]

f8 <- pmcmc(
            pomp(ou2,dprior=function (params, log, ...) {
              f <- sum(dnorm(params,mean=coef(ou2),sd=1,log=TRUE))
              if (log) f else exp(f)
            }),
            Nmcmc=1000,Np=500,verbose=FALSE,
            proposal=mvn.rw.adaptive(rw.sd=c(alpha.2=0.01,alpha.3=0.01),
              size.start=100,shape.start=100))
f8 <- continue(f8,Nmcmc=500,verbose=FALSE)
plot(f8)
require(coda)
trace <- window(conv.rec(f8,c("alpha.2","alpha.3")),start=500)
rejectionRate(trace)
effectiveSize(trace)
autocorr.diag(trace)
trace <- window(trace,thin=50)
plot(trace)
heidel.diag(trace)
geweke.diag(trace)

dev.off()
