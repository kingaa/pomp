library(pomp)
png(filename="steps-%02d.png",res=100)

pompExample(ricker)

set.seed(54588699L)

coef(ricker,"sigma") <- 0
tm <- sort(runif(n=20,max=3))
x <- trajectory(ricker,times=tm)["N",,]
y <- simulate(ricker,times=tm,states=TRUE)["N",,]
stopifnot(identical(x,y))

pompExample(euler.sir)
tm <- sort(runif(n=100,max=1))
x <- trajectory(euler.sir,times=tm)["I",,]
y <- simulate(euler.sir,times=tm,states=TRUE,seed=1232934371)["I",,]
table(cut(x-y,breaks=c(-Inf,seq(-0.2,0.2,by=0.01),Inf),ordered=T))

pompExample(ricker)
try(pomp(ricker,initializer=Csnippet("e=0; N = N_0; step = 0;")))
ricker <- pomp(ricker,
               rprocess=euler.sim(
                   Csnippet("double dW = rnorm(0,sqrt(dt));
                             N += r*N*(1-N)*dt+sigma*dW;
                             e += dW;
                             step += 1;"),
                   delta.t=0.1),
               skeleton=NULL,
               initializer=Csnippet("e=0; N = N_0; step = 0;"),
               zeronames="step",
               params=c(r=0.5,N_0=0.5,sigma=0.1,phi=10,c=1),
               paramnames=c("sigma","r","N_0"),statenames=c("N","e","step"))
time(ricker) <- sort(runif(n=20,max=50))
ricker <- simulate(ricker)
plot(ricker)

dev.off()
