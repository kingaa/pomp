## generate a bifurcation diagram for the Ricker map
p <- parmat(coef(ricker()),nrep=500)
p["r",] <- exp(seq(from=1.5,to=4,length=500))
x <- trajectory(ricker(),times=seq(from=1000,to=2000,by=1),params=p)
matplot(p["r",],x["N",,],pch='.',col='black',xlab="log(r)",ylab="N",log='x')
