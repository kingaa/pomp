library(pomp)

pompExample(bbs)

set.seed(48857734L)

coef(bbs)
coef(bbs,transform=TRUE)

bbs <- pomp(bbs,
            rprior=function(params,...){
              params["beta"] <- exp(runif(n=1,min=1,max=2))
              params["sigma"] <- runif(n=1,min=2,max=4)
              params
            }
            )

pf <- pfilter(bbs,Np=1000)

fit2 <- bsmc2(bbs,params=coef(bbs),Np=5000,transform=TRUE,
              est=c("beta","sigma"),smooth=0.2)
signif(coef(fit2),3)

fit3 <- traj.match(bbs,est=c("beta","sigma"),transform=TRUE)
signif(coef(fit3),3)

sim1 <- simulate(bbs,nsim=20,as.data.frame=TRUE,include.data=TRUE)
sim2 <- simulate(bbs,nsim=20,as.data.frame=TRUE,obs=TRUE,include.data=TRUE)
sim3 <- simulate(bbs,nsim=20,as.data.frame=TRUE,states=TRUE,include.data=TRUE)
