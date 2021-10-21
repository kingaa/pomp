options(digits=3)
png(filename="ou2-%02d.png",res=100)

library(pomp)

ou2() -> ou

set.seed(1438408329L)

plot(ou)
rinit(ou)
coef(ou)

stopifnot(all.equal(coef(ou),partrans(ou,coef(ou,transform=TRUE),dir="from")))
plot(s <- simulate(ou,seed=1438408329L))
pf <- freeze(pfilter(ou,Np=1000),seed=1438408329L)
plot(pf)
kf <- freeze(enkf(ou,Np=1000),seed=1438408329L)
plot(kf)
Kf <- kalmanFilter(
  ou,
  X0=rinit(ou),
  A=matrix(coef(ou,c("alpha_1","alpha_2","alpha_3","alpha_4")),2,2),
  Q={
    q <- matrix(c(coef(ou,c("sigma_1","sigma_2")),0,coef(ou,"sigma_3")),2,2)
    tcrossprod(q)
  },
  C=diag(2),
  R=diag(coef(ou,"tau")^2,2)
)
tj <- trajectory(ou,format="a")
matplot(time(ou),t(tj[,1,]),type="l",ylab="")

d <- dprocess(s,x=states(s),params=coef(s),times=time(s),log=TRUE)
plot(d[,],ylab="log prob")

try(dprocess(s,x=states(s)[,c(1:9,15)],params=coef(s),times=time(s)[c(1:9,15)]))

e <- emeasure(s,x=states(s),params=coef(s),times=time(s))
matplot(time(s),t(obs(s)),xlab="",ylab="",col=1:2)
matlines(time(s),t(apply(e,c(1,3),mean)),col=1:2)

v <- vmeasure(s,x=states(s),params=coef(s),times=time(s))
stopifnot(
  dim(v)==c(2,2,1,100),
  v[1,1,,]==v[2,2,,],
  v[1,2,,]==v[2,1,,],
  rownames(v)==c("y1","y2"),
  colnames(v)==rownames(v)
)

dev.off()
