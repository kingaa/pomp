library(pomp)

set.seed(398585L)
data(ou2)

time(ou2) <- 1:10

Np <- 10000

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
smc <- bsmc(
            ou2,
            params=prior,
            est=estnames,
            ntries=5,
            smooth=0.02,
            lower=prior.bounds[estnames,"lower"],
            upper=prior.bounds[estnames,"upper"]
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

smc <- bsmc(ou2,ntries=5,Np=5000,smooth=0.1,est=estnames,seed=455485L)
print(smc$eff.sample.size)
print(smc$log.evidence)
