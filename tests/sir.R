library(pomp)

png(filename='sir%02d.pdf',res=100)

tbasis <- seq(0,25,by=1/52)
basis <- periodic.bspline.basis(tbasis,nbasis=3,names="seas%d")

## some parameters
params <- c(
    gamma=26,mu=0.02,iota=0.01,
    beta1=400,beta2=480,beta3=320,
    beta.sd=1e-3,
    pop=2.1e6,
    rho=0.6,
    S.0=26/400,I.0=0.001,R.0=1-26/400
)

## set up the pomp object
## the C codes "sir_euler_simulator" and "sir_euler_density" are included in the "examples" directory (file "sir.c")
po <- pomp(
    times=1/52*seq.int(length=4*52),
    data=rbind(reports=numeric(52*4)),
    params=params,
    t0=0,
    tcovar=tbasis,
    covar=basis,
    zeronames=c("cases"),
    rprocess=euler.sim(
        delta.t=1/52/20,
        step.fun=function(t,x,params,covars,delta.t,...) {
            with(
                as.list(c(x,params)),
                {
                    beta <- sum(c(beta1,beta2,beta3)*covars)
                    dW <- rgammawn(n=1,sigma=beta.sd,dt=delta.t)
                    foi <- (iota+beta*I*dW/delta.t)/pop
                    trans <- c(
                        rpois(n=1,lambda=mu*pop*delta.t),
                        reulermultinom(n=1,size=S,rate=c(foi,mu),dt=delta.t),
                        reulermultinom(n=1,size=I,rate=c(gamma,mu),dt=delta.t),
                        reulermultinom(n=1,size=R,rate=c(mu),dt=delta.t)
                    )
                    c(
                        S=S+trans[1]-trans[2]-trans[3],
                        I=I+trans[2]-trans[4]-trans[5],
                        R=R+trans[4]-trans[6],
                        cases=cases+trans[4],
                        W=if (beta.sd>0) W+(dW-delta.t)/beta.sd else W,
                        B=trans[1],
                        SI=trans[2],
                        SD=trans[3],
                        IR=trans[4],
                        ID=trans[5],
                        RD=trans[6],
                        dW=dW
                    )
                }
            )
        }
    ),
    dprocess=onestep.dens(
        dens.fun=function(t1,t2,params,x1,x2,covars,...) {
            with(
                as.list(params),
                {
                    dt <- t2-t1
                    beta <- sum(c(beta1,beta2,beta3)*covars)
                    beta.var <- beta.sd^2
                    dW <- x2['dW']
                    foi <- (iota+beta*x1["I"]*dW/dt)/pop
                    probs <- c(
                        dpois(x=x2["B"],lambda=mu*pop*dt,log=T),
                        deulermultinom(x=x2[c("SI","SD")],size=x1["S"],rate=c(foi,mu),dt=dt,log=T),
                        deulermultinom(x=x2[c("IR","ID")],size=x1["I"],rate=c(gamma,mu),dt=dt,log=T),
                        deulermultinom(x=x2["RD"],size=x1["R"],rate=c(mu),dt=dt,log=T),
                        dgamma(x=dW,shape=dt/beta.var,scale=beta.var,log=T)
                    )
                    sum(probs)
                }
            )
        }
    ),
    skeleton=vectorfield(
        function(x,t,params,covars,...) {
            xdot <- rep(0,length(x))
            with(
                as.list(c(x,params)),
                {
                    beta <- sum(c(beta1,beta2,beta3)*covars)
                    foi <- (iota+beta*I)/pop
                    terms <- c(
                        mu*pop,
                        foi*S,
                        mu*S,
                        gamma*I,
                        mu*I,
                        mu*R
                    )
                    xdot[1:4] <- c(
                        terms[1]-terms[2]-terms[3],
                        terms[2]-terms[4]-terms[5],
                        terms[4]-terms[6],
                        terms[4]
                    )
                    xdot
                }
            )
        }),
    rmeasure=Csnippet("
        reports = nearbyint(rnorm(rho*cases,sqrt(rho*(1-rho)*cases)));
        if (reports < 0) reports = 0;
        "),
    paramnames="rho",statenames="cases",
    dmeasure=function(y,x,t,params,log,covars,...){
        with(
            as.list(c(x,params)),
            {
                if (y > 0) 
                    f <- diff(pnorm(q=y+c(-0.5,0.5),mean=rho*cases,sd=sqrt(rho*(1-rho)*cases),lower.tail=TRUE,log.p=FALSE))
                else
                    f <- pnorm(q=0.5,mean=rho*cases,sd=sqrt(rho*(1-rho)*cases),lower.tail=TRUE,log.p=FALSE)
                if (log) log(f) else f
            }
        )
    },
    initializer=function(params,t0,...){
        with(
            as.list(params),
            {
                fracs <- c(S.0,I.0,R.0)
                x0 <- c(
                    round(pop*fracs/sum(fracs)), # make sure the three compartments sum to 'pop' initially
                    rep(0,9)	# zeros for 'cases', 'W', and the transition numbers
                )
                names(x0) <- c("S","I","R","cases","W","B","SI","SD","IR","ID","RD","dW")
                x0["cases"] <- gamma/52*x0["I"]
                x0
            }
        )
    }
)

show(po)

x <- pomp(po,measurement.model=reports~binom(size=cases,prob=rho),
          rmeasure=function(x,t,params,covars,...){1})

set.seed(3049953)
## simulate from the model
x <- simulate(po,nsim=3)

t1 <- seq(0,4/52,by=1/52/25)
X1 <- simulate(po,nsim=10,states=TRUE,obs=TRUE,times=t1)

t2 <- seq(0,2,by=1/52)
X2 <- simulate(po,nsim=1,states=TRUE,obs=TRUE,times=t2)

t3 <- seq(0,20,by=1/52)
X3 <- trajectory(po,times=t3,hmax=1/52)

f1 <- dprocess(
    po,
    x=X1$states[,,31:40],
    times=t1[31:40],
    params=parmat(params,nrep=10),
    log=TRUE
)
print(apply(f1,1,sum),digits=4)

g1 <- dmeasure(
    po,
    y=rbind(reports=X1$obs[,7,]),
    x=X1$states,
    times=t1,
    params=parmat(params,nrep=10),
    log=TRUE
)
print(apply(g1,1,sum),digits=4)

h1 <- skeleton(
    po,
    x=X2$states[,1,55:70,drop=FALSE],
    t=t2[55:70],
    params=params
)
print(h1[c("S","I","R"),,],digits=4)

## now repeat using the compiled native codes built into the package
pompExample(euler.sir)
po <- euler.sir

stopifnot(all.equal(
    partrans(po,coef(po,trans=FALSE),dir="to"),
    coef(po,trans=TRUE)))
stopifnot(all.equal(
    partrans(po,coef(po,trans=TRUE),dir="from"),
    coef(po,trans=FALSE),
    tolerance=1.5e-6))

set.seed(3049953)
## simulate from the model
x <- simulate(po,nsim=100)
t3 <- seq(0,20,by=1/52)
X4 <- trajectory(po,times=t3,hmax=1/52)

g2 <- dmeasure(
    po,
    y=rbind(reports=X1$obs[,7,]),
    x=X1$states,
    times=t1,
    params=parmat(coef(po),nrep=10),
    log=TRUE
)
print(apply(g2,1,sum),digits=4)

h2 <- skeleton(
    po,
    x=X2$states[,1,55:70,drop=FALSE],
    t=t2[55:70],
    params=coef(po)
)
print(h2[c("S","I","R"),,],digits=4)

stopifnot(max(abs(g2-g1),na.rm=T)<.Machine$double.eps*100)
stopifnot(max(abs(h2-h1),na.rm=T)<.Machine$double.eps*2^16)

invisible(states(po)[,1:2])
time(po) <- seq(0,1,by=1/52)
invisible(states(po)[,1:3])
invisible(states(simulate(po))[,1:3])

po <- window(euler.sir,start=1,end=2)
timezero(po)
timezero(po)<-2*time(po)[1]-time(po)[2]

dev.off()

## test of vectorfield integrator

pompExample(euler.sir)

po <- pomp(
    window(euler.sir,end=2),
    skeleton=vectorfield(
        function(x,t,params,...) {
            xdot <- rep(0,length(x))
            with(
                as.list(c(x,params)),
                {
                    covars <- as.numeric(periodic.bspline.basis(t,nbasis=3,degree=3,period=1))
                    beta <- sum(c(beta1,beta2,beta3)*covars)
                    foi <- (iota+beta*I)/pop
                    terms <- c(
                        mu*pop,
                        foi*S,
                        mu*S,
                        gamma*I,
                        mu*I,
                        mu*R
                    )
                    xdot[1:4] <- c(
                        terms[1]-terms[2]-terms[3],
                        terms[2]-terms[4]-terms[5],
                        terms[4]-terms[6],
                        terms[4]
                    )
                    xdot
                }
            )
        }
    )
)

x1 <- trajectory(po,hmax=1/52,as.data.frame=T)
x2 <- trajectory(window(euler.sir,end=2),hmax=1/52,as.data.frame=T)

stopifnot(identical(round(x1$S,3),round(x2$S,3)))
stopifnot(identical(round(x1$I,3),round(x2$I,3)))
stopifnot(identical(round(x1$cases,3),round(x2$cases,3)))

pomp(euler.sir,
     rmeasure=Csnippet("
       const SEXP bob = get_pomp_userdata(\"bob\");
       if (*INTEGER(bob) != 33) {
         Rprintf(\"error!\");
       }
       double mean = rho*cases;
       double sd = sqrt(cases*rho*(1-rho));
       reports = nearbyint(rnorm(mean,sd));
       reports = (reports > 0) ? reports : 0;"),
    statenames=c("cases"),
    paramnames=c("rho"),bob=33L) -> po

simulate(po) -> po
