library(pomp)

data(ricker)

pdf(file="ricker.pdf")

tj.1 <- trajectory(ricker)
plot(time(ricker),tj.1[1,,],type='l')
tj.2 <- trajectory(ricker,times=c(30:50),t0=0)
lines(30:50,tj.2[1,,],col='red',lwd=2)
max(abs(tj.1[,,time(ricker)>=30]-tj.2[,,]))

tj.3 <- trajectory(ricker,as.data.frame=TRUE)
plot(tj.3)
tj.3 <- trajectory(ricker,as.data.frame=TRUE,params=parmat(coef(ricker),3),times=1:100)
plot(N~time,data=tj.3,subset=traj==3,type='l')

sm <- simulate(ricker,seed=343995,as.data.frame=TRUE)
sm1 <- as.data.frame(simulate(ricker,seed=343995))
stopifnot(max(abs(as.matrix(sm[names(sm1)])-as.matrix(sm1)))==0)

sm <- simulate(ricker,nsim=3,seed=343995,as.data.frame=TRUE)
print(names(sm))
print(dim(sm))

sm1 <- simulate(ricker,nsim=3,obs=T,seed=343995,as.data.frame=TRUE)
print(names(sm1))
print(dim(sm1))
stopifnot(max(abs(as.matrix(sm[names(sm1)])-as.matrix(sm1)))==0)

sm1 <- simulate(ricker,nsim=3,states=T,seed=343995,as.data.frame=TRUE)
print(names(sm1))
print(dim(sm1))
stopifnot(max(abs(as.matrix(sm[names(sm1)])-as.matrix(sm1)))==0)

sm1 <- simulate(ricker,nsim=3,states=T,obs=T,seed=343995,as.data.frame=TRUE)
print(names(sm1))
print(dim(sm1))
stopifnot(max(abs(as.matrix(sm[names(sm1)])-as.matrix(sm1)))==0)

sm <- simulate(ricker,nsim=1,states=T,obs=T,seed=343995,as.data.frame=TRUE)
sm1 <- as.data.frame(simulate(ricker,seed=343995))
print(names(sm))
print(dim(sm))
stopifnot(max(abs(as.matrix(sm[names(sm1)])-as.matrix(sm1)))==0)

po <- ricker
try(
    coef(po,"r")
    )
coef(po,c("r","phi")) <- c(0,0)
coef(po,c("log.r","phi")) <- c(a=0,b=1)
coef(po,c("log.r","phi")) <- 1
coef(po) <- c(phi=1,log.r=3.5,N.0=10,e.0=0,sigma=0)
coef(po)
coef(po,"new") <- 3
plot(simulate(po))
coef(po)

dev.off()
