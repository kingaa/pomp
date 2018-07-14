library(pomp)

pompExample(bbs)

set.seed(48857734L)

coef(bbs)
coef(bbs,transform=TRUE)

bbs <- pomp(bbs,rprior=Csnippet("
              Beta = exp(runif(1,2));
              sigma = runif(2,4);"),
            paramnames=c("Beta","sigma"))

pf <- pfilter(bbs,Np=1000)

fit2 <- bsmc2(bbs,params=coef(bbs),Np=5000,transform=TRUE,
              est=c("Beta","sigma"),smooth=0.2)
signif(coef(fit2),3)

fit3 <- traj.match(bbs,est=c("Beta","sigma"),transform=TRUE)
signif(coef(fit3),3)

sim1 <- simulate(bbs,nsim=20,as.data.frame=TRUE,include.data=TRUE)
sim2 <- simulate(bbs,nsim=20,as.data.frame=TRUE,obs=TRUE,include.data=TRUE)
sim3 <- simulate(bbs,nsim=20,as.data.frame=TRUE,states=TRUE,include.data=TRUE)
