require(pomp)

flu <- read.csv2(text="
day;reports
1;3
2;8
3;28
4;76
5;222
6;293
7;257
8;237
9;192
10;126
11;70
12;28
13;12
14;5
")

po <- pomp(
           data=flu,
           times="day",
           t0=0,
           rprocess=euler.sim(
             step.fun="_sir_euler_simulator",
             delta.t=1/12,
             PACKAGE="pomp"
             ),
           skeleton.type="vectorfield",
           skeleton="_sir_ODE",
           measurement.model=reports~norm(mean=rho*cases,sd=1+sigma*cases),
           PACKAGE="pomp",
           obsnames = c("reports"),
           statenames=c("S","I","R","cases","W"),
           paramnames=c(
             "gamma","mu","iota",
             "beta","beta.sd","pop","rho",
             "S.0","I.0","R.0"
             ),
           zeronames=c("cases"),
           comp.names=c("S","I","R"),
           ic.names=c("S.0","I.0","R.0"),
           nbasis=1L,
           degree=0L,
           period=1.0,
           logvar=c(
             "beta","gamma","mu","iota","sigma","beta.sd",
             "S.0","I.0","R.0"
             ),
           logitvar="rho",
           parameter.inv.transform=function (params, logvar, logitvar, ...) {
             params[logvar] <- log(params[logvar])
             params[logitvar] <- qlogis(params[logitvar])
             params
           },
           parameter.transform=function (params, logvar, logitvar, ...) {
             params[logvar] <- exp(params[logvar])
             params[logitvar] <- plogis(params[logitvar])
             params
           },
           initializer=function(params, t0, comp.names, ic.names, ...) {
             snames <- c("S","I","R","cases","W")
             fracs <- params[ic.names]
             x0 <- numeric(length(snames))
             names(x0) <- snames
             x0[comp.names] <- round(params['pop']*fracs/sum(fracs))
             x0
           }
           )

coef(po) <- c(
              gamma=1/3,mu=0.0,iota=0.0,
              beta=1.4,
              beta.sd=0,
              pop=1400,
              rho=0.9,sigma=3.6,
              S.0=0.999,I.0=0.001,R.0=0
              )

bbs <- po

assign("bbs",bbs,envir=.GlobalEnv)
c("bbs")
