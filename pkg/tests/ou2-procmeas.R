library(pomp)

data(ou2)

po <- window(ou2,end=10)

set.seed(3434388L)

pmat <- do.call(cbind,rep(list(coef(po)),3))
sims <- simulate(po,states=T,obs=T,params=pmat)

dp <- dprocess(po,x=sims$states,times=time(po),params=pmat,log=T)
dm <- dmeasure(po,x=sims$states,y=obs(po),times=time(po),params=pmat,log=T)

apply(dp,1,sum)
apply(dm,1,sum)
