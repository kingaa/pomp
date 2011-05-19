library(pomp)

data(ou2)

set.seed(3434388L)

pmat <- do.call(cbind,rep(list(coef(ou2)),3))
sims <- simulate(ou2,states=T,obs=T,params=pmat)

dp <- dprocess(ou2,x=sims$states,times=time(ou2),params=pmat,log=T)
dm <- dmeasure(ou2,x=sims$states,y=obs(ou2),times=time(ou2),params=pmat,log=T)

apply(dp,1,sum)
apply(dm,1,sum)
