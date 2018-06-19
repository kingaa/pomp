library(pomp)

set.seed(93885485L)

pompExample(ou2)
true.p <- coef(ou2)
simdata <- simulate(ou2,nsim=5,params=true.p,seed=394885)
guess.p <- true.p

x <- ou2
coef(x) <- numeric()
try(trajectory(x))
try(trajectory(pomp(ou2,skeleton=NULL)))
pompExample(euler.sir)
try(trajectory(euler.sir,atol=-100))

x <- sapply(
    simdata,
    function (d) {
        res <- traj.match(
            d,
            start=guess.p,
            est=c('alpha.1','alpha.4','x1.0','x2.0','tau'),
            maxit=2000,
            reltol=1e-8
        )
        c(conv=res$convergence,loglik=logLik(res),coef(res))
    }
)
range(x['conv',])
range(x['loglik',])
print(
    cbind(
        truth=true.p,
        mean.est=apply(x[names(true.p),],1,mean),
        bias=apply(x[names(true.p),],1,mean)-true.p,
        se=apply(x[names(true.p),],1,sd)
    ),
    digits=4
)

summary(traj.match(ou2,est=c('alpha.1','alpha.4','x1.0','x2.0','tau'),method="subplex",maxit=100))

invisible(capture.output(summary(traj.match(ou2,est=c('alpha.1','x1.0','alpha.4','x2.0','tau'),method="sannbox",trace=4,parscale=0.1,maxit=100))))

try(traj.match(ou2,est=c('alpha.1','x1.0','alpha.4','x2.0','tau'),method="sannbox",parscale=0.1,maxit=100,sched=rep(1,10)))

summary(traj.match(ou2))
ofun <- traj.match.objfun(ou2)
try(traj.match.objfun(ou2,params=c("A")))
ofun <- traj.match.objfun(ou2,est=c('x1.0','x2.0','alpha.1','alpha.4','tau'))
print(ofun(coef(ou2,c('x1.0','x2.0','alpha.1','alpha.4','tau'))))

try(ofun <- traj.match.objfun(ou2,est=c('x1.0','x2.0','bob','alpha.4','harry')))
try(traj.match(ou2,est=c('x1.0','x2.0','bob','alpha.4','tau')))

fit <- optim(par=c(0,0,1,0.5,0.5),fn=ofun,method="Nelder-Mead",control=list(maxit=400))
print(fit)
stopifnot(fit$convergence==0)

pompExample(ou2)
p <- coef(ou2)
ou2@params <- numeric(0)
res <- traj.match(ou2,start=p)
print(coef(res),digits=4)
res <- traj.match(res,est=c('alpha.1','alpha.4','x1.0','x2.0','tau'),maxit=2000,reltol=1e-8)
print(coef(res),digits=4)
print(p,digits=4)
traj.match(res) -> res

pomp(ou2,skeleton=map(function (x, t, params, ...) {unname(x+1)})) -> po
trajectory(po,params=p,t0=-10,as.data.frame=TRUE) -> x
sapply(x[,1:2],range)

pp <- p
pp[c("alpha.1","alpha.2","x1.0","x2.0","tau")] <- c(1,-1,-5,5,2)
traj.match(ou2,start=pp,est=c("alpha.1","alpha.2","x1.0","x2.0","tau"),
           method="nloptr",algorithm="NLOPT_GN_MLSL",
           lb=c(alpha.1=0,alpha.2=-1,x1.0=-10,x2.0=0,tau=0),
           ub=c(alpha.1=2,alpha.2=0,x1.0=0,x2.0=10,tau=5),
           xtol_rel=1e-6,maxeval=1000,
           local_opts=list(algorithm="NLOPT_LN_NELDERMEAD",
                           xtol_rel=1e-3)) -> tm
stopifnot(summary(tm)$convergence==5)
traj.match(ou2,start=pp,est=c("alpha.1","alpha.2","x1.0","x2.0","tau"),
           method="nloptr",algorithm="NLOPT_GN_MLSL",
           lower=c(0,-1,-10,0,0),
           upper=c(2,0,0,10,5),
           xtol_rel=1e-6,maxeval=1000,
           local_opts=list(algorithm="NLOPT_LN_NELDERMEAD",
                           xtol_rel=1e-3)) -> tm
stopifnot(summary(tm)$convergence==5)

pompExample(ricker)
try(traj.match(ricker,start=numeric(0)))
try(traj.match(ricker,transform=TRUE,est=c("bob")))
