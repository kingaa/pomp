library(pomp)

pompExample("ricker")

pdf(file="ricker.pdf")

tj.1 <- trajectory(ricker)
plot(time(ricker),tj.1[1,,],type='l')
tj.2 <- trajectory(ricker,times=c(30:50),t0=0)
lines(30:50,tj.2[1,,],col='red',lwd=2)
stopifnot(max(abs(tj.1[,,time(ricker)>=30]-tj.2[,,]))==0)

tj.3 <- trajectory(ricker,as.data.frame=TRUE)
plot(tj.3)
tj.3 <- trajectory(ricker,as.data.frame=TRUE,params=parmat(coef(ricker),3),times=1:100)
plot(N~time,data=tj.3,subset=traj==3,type='l')

sm <- as.data.frame(simulate(ricker,seed=343995))
sm1 <- simulate(ricker,seed=343995,as.data.frame=TRUE)
stopifnot(max(abs(as.matrix(sm1[names(sm)])-as.matrix(sm)))==0)
sm1 <- simulate(ricker,seed=343995,states=TRUE,obs=TRUE,as.data.frame=TRUE)
stopifnot(max(abs(as.matrix(sm1[names(sm)])-as.matrix(sm)))==0)

sm1 <- simulate(ricker,nsim=3,seed=343995,as.data.frame=TRUE)
stopifnot(all(names(sm1)==c("time","y","N","e","sim")))
stopifnot(all(dim(sm1)==c(153,5)))

sm1 <- simulate(ricker,nsim=3,states=T,seed=343995,as.data.frame=TRUE)
stopifnot(all(names(sm1)==c("N","e","sim","time")))
stopifnot(all(dim(sm1)==c(153,4)))

sm1 <- simulate(ricker,nsim=3,obs=T,seed=343995,as.data.frame=TRUE)
stopifnot(all(names(sm1)==c("y","sim","time")))
stopifnot(all(dim(sm1)==c(153,3)))

po <- ricker
try(
    coef(po,"log.r")
    )
coef(po,c("r","phi")) <- c(0,0)
coef(po,c("r","phi")) <- c(a=0,b=1)
coef(po,c("r","phi")) <- 1
coef(po) <- c(phi=1,r=3.5,N.0=10,e.0=0,sigma=0)
coef(po)
coef(po,"new") <- 3
plot(simulate(po))
coef(po)

dev.off()
