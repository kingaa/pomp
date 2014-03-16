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

fit1 <- bsmc(bbs,params=coef(bbs),Np=1000,transform=TRUE,est=c("beta","sigma"),smooth=0.2)
signif(coef(fit1),3)

fit2 <- traj.match(bbs,est=c("beta","sigma"),transform=TRUE)
signif(coef(fit2),3)
