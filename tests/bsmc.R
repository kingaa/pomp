library(pomp)

options(digits=2)

png(filename="bsmc-%02d.png",res=100)

set.seed(398585L)
pompExample(ou2)

time(ou2) <- 1:10

Np <- 10000

try(smc <- bsmc(ou2,Np=2,smooth=0.01,est=estnames,
                 tol=1e-2,max.fail=100))

prior.bounds <- rbind(
                      alpha.2=c(-0.55,-0.45),
                      alpha.3=c(0.25,0.35)
                      )
colnames(prior.bounds) <- c("lower","upper")

estnames <- rownames(prior.bounds)

prior <- matrix(data=coef(ou2),nrow=length(coef(ou2)),ncol=Np)
rownames(prior) <- names(coef(ou2))
try(bsmc(ou2,params=prior))
for (n in estnames) {
  prior[n,] <- runif(n=Np,min=prior.bounds[n,1],max=prior.bounds[n,2])
}
try(bsmc(ou2,est=c("alpha.2","alpha.3"),params={p <- prior; rownames(p) <- NULL; p}))
try(bsmc(ou2,params=prior,smooth=2))

try(bsmc(ou2,est=c("alpha.2","alpha.3"),Np=1,smooth=1e-100))

garb <- ou2
coef(garb) <- numeric(0)
try(garb <- bsmc(garb))

##Run Liu & West particle filter
bsmc(ou2,params=prior,smooth=0.02,seed=49959,Np=100) -> smc
smc <- bsmc(ou2,est="alpha.2",params=prior,smooth=0.02)
prior <- smc$prior
post <- smc$post

try(bsmc(ou2,params=prior,est=estnames,ntries=5,smooth=0.02,lower=0,upper=c(0,1)))
try(bsmc(ou2,params=prior,est=estnames,ntries=5,smooth=0.02,lower=-100,upper=c(111,33,222)))

smc <- bsmc(ou2,params=prior,est=estnames,ntries=5,smooth=0.02,
            lower=prior.bounds[estnames,"lower"],
            upper=prior.bounds[estnames,"upper"]
            )
prior <- smc$prior
post <- smc$post
print(min(smc$eff.sample.size))
print(smc$log.evidence)

ou2 <- pomp(ou2,
            rprior=function(params,...){
              params
            }
            )

smc <- bsmc(ou2,ntries=5,Np=5000,smooth=0.1,est=estnames)
smc <- bsmc(ou2,ntries=5,transform=TRUE,Np=5000,smooth=0.1,est=estnames)
plot(smc,y=NA)
try(plot(smc,pars=c("george","gracie")))
print(smc$eff.sample.size)
print(smc$log.evidence)

smc <- bsmc(ou2,ntries=1,transform=TRUE,Np=2,smooth=0.01,est=estnames,
            tol=1e-2,max.fail=100)

try(bsmc(pomp(ou2,dmeasure=function(y,x,t,params,log,...) stop("oof!")),
         Np=100,est=c("alpha.1","alpha.2"),transform=TRUE,smooth=0.2))
neval <- 0
try(bsmc(pomp(ou2,
              dmeasure=function(y,x,t,params,log,...) {
                  neval <<- neval+1
                  if (neval>100) stop("oof!") else if (log) 0 else 1
              }),
         Np=100,est=c("alpha.1","alpha.2"),transform=TRUE,smooth=0.2))

try(bsmc(ou2,Np=100,est=c("alpha.1","bob"),transform=TRUE,smooth=1))

dev.off()
