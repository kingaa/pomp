library(pomp)
set.seed(1066L)

pdf(file="probe.pdf")

pompExample(ou2)

pm.ou2 <- probe(
    ou2,
    probes=list(
        y1.mean=probe.mean(var="y1"),
        y2.mean=probe.mean(var="y2"),
        y1.sd=probe.sd(var="y1"),
        y2.sd=probe.sd(var="y2")
    ),
    nsim=500
)

pm.po <- probe(
    ou2,
    params=c(
        alpha.1=0.1,alpha.2=-0.5,alpha.3=0.3,alpha.4=0.2,
        sigma.1=3,sigma.2=-0.5,sigma.3=2,
        tau=1,
        x1.0=0,x2.0=0
    ),
    probes=list(
        y1.mean=probe.mean(var="y1"),
        y2.mean=probe.mean(var="y2"),
        y1.sd=probe.sd(var="y1"),
        y2.sd=probe.sd(var="y2")
    ),
    nsim=500
)

invisible(summary(pm.ou2))
invisible(summary(pm.po))

plot(pm.ou2)
plot(pm.po)

pm.ou2 <- probe(
    ou2,
    probes=list(
        y1acf=probe.acf(var="y1",lags=c(0,1,2),type="corr"),
        y2acf=probe.acf(var=c("y2"),lags=c(0,1,2)),
        y12ccf=probe.ccf(var=c("y2","y1"),lags=c(3,8))
    ),
    nsim=500
)
plot(pm.ou2)

pb <- probe(
    ou2,
    probes=list(
        y1=probe.quantile(var="y1",prob=seq(0.1,0.9,by=0.1)),
        probe.acf(var=c("y1","y2"),lags=c(0,1,4,7),transform=identity),
        pd=probe.period(var="y1",kernel.width=3)
    ),
    nsim=200
)
plot(pb)

po <- ou2
coef(po,c("alpha.2","alpha.3")) <- c(0,0)
coef(po,c("sigma.2","sigma.1","sigma.3")) <- c(0,0.0,0.0)
coef(po,c("tau")) <- c(0.0)
po <- simulate(po)
pb <- probe(
    po,
    probes=list(
        probe.acf(var=c("y1","y2"),lags=c(0,1),type="cor"),
        probe.nlar("y1",lags=1,powers=1),
        probe.nlar("y2",lags=1,powers=1)
    ),
    nsim=1000,
    seed=1066L
)
x <- as.data.frame(po)
mx <- sapply(x,mean)
x <- sweep(x,2,mx)
y1 <- head(x$y1,-1)
z1 <- tail(x$y1,-1)
y2 <- head(x$y2,-1)
z2 <- tail(x$y2,-1)
small.diff <- pb$datvals-c(mean(y1*z1)/mean(x$y1^2),mean(y2*z2)/mean(x$y2^2),mean(y1*z1)/mean(y1*y1),mean(y2*z2)/mean(y2*y2))
stopifnot(max(abs(small.diff))<.Machine$double.eps*100)

po <- simulate(ou2)
pb <- probe(
    po,
    probes=list(
        probe.acf(var=c("y1"),lags=c(0,1,2),type="cov"),
        probe.ccf(vars=c("y1","y1"),lags=c(0,1,2),type="cov")
    ),
    nsim=1000,
    seed=1066L
)
plot(pb)

pb <- probe(
    po,
    probes=probe.ccf(vars=c("y1","y2"),lags=c(-5,-3,1,4,8)),
    nsim=1000,
    seed=1066L
)
plot(pb)

pb <- probe(
    po,
    probes=probe.ccf(vars=c("y1","y2"),lags=c(-5,-3,1,4,8),type="corr"),
    nsim=1000,
    seed=1066L
)
plot(pb)

head(as(pb,"data.frame"))

pompExample(ou2)

good <- probe(
    ou2,
    probes=list(
        y1.mean=probe.mean(var="y1"),
        y2.mean=probe.mean(var="y2"),
        y1.sd=probe.sd(var="y1"),
        y2.sd=probe.sd(var="y2"),
        extra=function(x)max(x["y1",])
    ),
    nsim=500
)

ofun <- probe.match.objfun(ou2,est=c("alpha.1","alpha.2"),
                           probes=good$probes,nsim=100,
                           seed=349956868L
                           )

library(nloptr)
fit1 <- nloptr(
    coef(good,c("alpha.1","alpha.2")),
    ofun,
    opts=list(
        algorithm="NLOPT_LN_SBPLX",
        xtol_rel=1e-10,
        maxeval=1000
    )
)
fit2 <- probe.match(
    good,
    est=c("alpha.1","alpha.2"),
    nsim=100,
    algorithm="NLOPT_LN_SBPLX",
    xtol_rel=1e-10,
    maxeval=1000,
    seed=349956868L
)

all.equal(fit1$solution,unname(coef(fit2,fit2$est)))

pompExample(ricker)

set.seed(6457673L)
z <- as.numeric(obs(ricker))

po <- ricker

pb <- probe(
    po,
    probes=probe.median("y"),
    nsim=1000,
    seed=838775L
)
plot(pb)
invisible(summary(pb))

pb <- probe(
    po,
    probes=probe.nlar("y",lags=c(1,2,3),powers=c(1,1,1),transform="sqrt"),
    nsim=1000,
    seed=838775L
)
plot(pb)
invisible(summary(pb))

pb <- probe(
    po,
    probes=probe.nlar("y",lags=c(1,2,3),powers=1,transform="sqrt"),
    nsim=1000,
    seed=838775L
)
plot(pb)
invisible(summary(pb))

pb <- probe(
    po,
    probes=probe.nlar("y",lags=1,powers=c(1,2,3),transform="sqrt"),
    nsim=1000,
    seed=838775L
)
plot(pb)
invisible(summary(pb))

pb <- probe(
    po,
    probes=probe.marginal(
        var="y",
        transform=sqrt,
        ref=z,
        diff=1,
        order=3
    ),
    nsim=1000,
    seed=838775L
)
invisible(pb$datvals)
invisible(summary(pb))
plot(pb)

pb <- probe(
    po,
    probes=list(
        probe.marginal(
            var="y",
            transform=sqrt,
            ref=z,
            diff=1,
            order=3
        ),
        probe.acf(
            var="y",
            lags=c(0,1,3,5)
        ),
        mean=probe.mean(var="y",transform=sqrt)
    ),
    nsim=1000,
    seed=838775L
)
invisible(pb$datvals)
invisible(summary(pb))
plot(pb)

pbm <- probe.match(pb)
plot(pbm)
invisible(summary(pbm))

coef(po) <- c(r=10,sigma=0.3,phi=20,N.0=5,e.0=0)

pb <- probe(
    po,
    probes=probe.marginal(
        var="y",
        transform=sqrt,
        ref=z,
        diff=1,
        order=3
    ),
    nsim=1000,
    seed=838775L
)
invisible(pb$datvals)
invisible(summary(pb))
plot(pb)

pm <- probe.match(
    pb,
    est=c("r","phi","N.0"),
    transform=TRUE,
    parscale=c(0.1,0.1,0.1),
    nsim=1000,
    seed=838775L,
    method="Nelder-Mead",
    reltol=1e-7,
    fail.value=1e9
)
plot(pm)

invisible(cbind(truth=coef(ricker),est=coef(pm),guess=coef(po)))

pb <- probe(
    po,
    probes=probe.nlar(
        var="y",
        transform=sqrt,
        lags=1,
        powers=c(1,2,3)
    ),
    nsim=1000,
    seed=838775L
)
invisible(pb$datvals)
invisible(summary(pb))
plot(pb)

pm <- probe.match(
    pb,
    est=c("r","phi","N.0"),
    transform=TRUE,
    parscale=c(0.1,0.1,0.1),
    nsim=1000,
    seed=838775L,
    method="Nelder-Mead",
    reltol=1e-7,
    fail.value=1e9
)
plot(pm)
plot(as(pm,"pomp"),variables="y")
pm <- probe.match(pm,seed=613980375L)
plot(simulate(pm),variables="y")

invisible(cbind(truth=coef(ricker),est=coef(pm),guess=coef(po)))

pb <- probe(
    po,
    probes=probe.marginal(
        var="y",
        transform=sqrt,
        ref=runif(length(time(ricker))),
        diff=2,
        order=3
    ),
    nsim=1000,
    seed=838775L
)
invisible(pb$datvals)
invisible(summary(pb))
plot(pb)

pb <- probe(
    ricker,
    probes=probe.acf(
        var="y",
        transform=sqrt,
        lags=seq.int(from=0,to=5),
        type="cov"
    ),
    nsim=1000,
    seed=838775L
)
invisible(pb$datvals)
invisible(summary(pb))
plot(pb)

pb <- probe(
    ricker,
    probes=probe.acf(
        var="y",
        transform=sqrt,
        lags=seq.int(from=1,to=5),
        type="cor"
    ),
    nsim=1000,
    seed=838775L
)
invisible(pb$datvals)
invisible(summary(pb))
plot(pb)

pb <- probe(
    ricker,
    probes=list(
        v=probe.var("y",transform=sqrt),
        probe.acf(
            var="y",
            transform=sqrt,
            lags=c(0,1,2),
            type="cov"
        ),
        probe.nlar(
            var="y",
            transform=sqrt,
            lags=c(1,2),
            powers=1
        )
    ),
    nsim=1000,
    seed=838775L
)
invisible(pb$datvals)
invisible(summary(pb))
plot(pb)

try(
    probe(
        ricker,
        probes=list(
            mn=probe.mean("y",transform=sqrt,trim=0.1),
            md=probe.median("y",na.rm=FALSE),
            wacko=function(y)y[sample.int(n=length(y),size=sample.int(n=3,size=1))]
        ),
        nsim=100,
        seed=838775L
    )
)

try(
    probe(
        ricker,
        probes=list(
            mn=probe.mean("y",transform=sqrt,trim=0.1),
            md=function(y)median(as.numeric(y)),
            wacko=function(y) if (y[1]==68) 1 else c(1,2)
        ),
        nsim=100,
        seed=838775L
    )
)


try(
    probe(
        ricker,
        probes=list(
            mn=probe.mean("y",transform=sqrt,trim=0.1),
            md=function(y)median(as.numeric(y)),
            wacko=function(y) if (y[28]==107) c(1,2) else 1
        ),
        nsim=20,
        seed=838775L
    )
)

dev.off()
