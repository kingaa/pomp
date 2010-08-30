library(pomp)

data(ou2)

## fix some parameters
p <- c(
       alpha.1=0.9,alpha.2=0,alpha.3=0,alpha.4=0.99,
       sigma.1=1,sigma.2=0,sigma.3=2,
       tau=1,
       x1.0=50,x2.0=-50       
       )

tic <- Sys.time()
ou2.sim <- simulate(ou2,params=p,nsim=100,seed=32043858)
toc <- Sys.time()
print(toc-tic)

coef(ou2,c('x1.0','x2.0')) <- c(-50,50)

ou2.sim <- simulate(ou2)
x <- simulate(ou2,nsim=3,states=T)
y <- simulate(ou2,nsim=3,obs=T)
z <- simulate(ou2,nsim=3,obs=T,states=T)
