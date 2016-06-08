library(pomp)

pdf(file="bsmc2.pdf")

set.seed(398585L)
pompExample(ou2)

time(ou2) <- 1:10

Np <- 50000

prior.bounds <- rbind(
                      alpha.2=c(-0.55,-0.45),
                      alpha.3=c(0.25,0.35)
                      )
colnames(prior.bounds) <- c("lower","upper")

estnames <- rownames(prior.bounds)

prior <- matrix(data=coef(ou2),nrow=length(coef(ou2)),ncol=Np)
rownames(prior) <- names(coef(ou2))
for (n in estnames) {
  prior[n,] <- runif(n=Np,min=prior.bounds[n,1],max=prior.bounds[n,2])
}

##Run Liu & West particle filter
tic <- Sys.time()
smc <- bsmc2(
             ou2,
             params=prior,
             est=estnames,
             smooth=0.02
             )
toc <- Sys.time()

prior <- smc$prior
post <- smc$post

print(etime <- toc-tic)

print(
      cbind(
            prior.mean=apply(prior,1,mean),
            posterior.mean=apply(post,1,mean),
            truth=coef(ou2),
            t(apply(post,1,quantile,c(0.025,0.5,0.975)))
            )
      )

print(min(smc$eff.sample.size))
print(smc$log.evidence)

ou2 <- pomp(ou2,
            rprior=function(params,...){
              params
            }
            )

smc <- bsmc2(ou2,Np=25000,smooth=0.1,est=estnames)
print(smc$eff.sample.size)
print(smc$log.evidence)

pompExample(ricker)

set.seed(6457673L)

po <- pomp(
           ricker,
           rprior=function (params, ...) {
             params["r"] <- exp(runif(n=1,min=2,max=5))
             params["sigma"] <- runif(n=1,min=0.1,max=1)
             params
           }
           )

Np <- 10000

fit <- bsmc2(po,Np=100,est=c("r","sigma"),transform=TRUE,smooth=0.2)

invisible(apply(fit$prior[c("r","sigma"),],1,mean))

invisible(apply(fit$post[c("r","sigma"),],1,mean))

invisible(coef(fit))

plot(fit,thin=300)

dev.off()
