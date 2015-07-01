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

pomp(
     data=flu,
     times="day",
     t0=0,
     params=c(
       gamma=1/3,mu=0.0,iota=0.0,
       beta=1.4,
       beta.sd=0,
       pop=1400,
       rho=0.9,sigma=3.6,
       S_0=0.999,I_0=0.001,R_0=0
       ),
     rprocess=euler.sim(
       step.fun="_sir_euler_simulator",
       delta.t=1/12,
       PACKAGE="pomp"
       ),
     skeleton.type="vectorfield",
     skeleton="_sir_ODE",
     measurement.model=reports~nbinom(mu=rho*cases,size=1/sigma^2),
     PACKAGE="pomp",
     statenames=c("S","I","R","cases","W"),
     paramnames=c(
       "gamma","mu","iota",
       "beta","beta.sd","pop","rho",
       "S_0","I_0","R_0"
       ),
     zeronames=c("cases"),
     nbasis=1L,
     degree=0L,
     period=1.0,
     logvar=c(
       "beta","gamma","mu","iota","sigma","beta.sd",
       "S_0","I_0","R_0"
       ),
     logitvar="rho",
     toEstimationScale=function (params, logvar, logitvar, ...) {
       params[logvar] <- log(params[logvar])
       params[logitvar] <- qlogis(params[logitvar])
       params
     },
     fromEstimationScale=function (params, logvar, logitvar, ...) {
       params[logvar] <- exp(params[logvar])
       params[logitvar] <- plogis(params[logitvar])
       params
     },
     initializer="_sir_init"
     ) -> bbs

c("bbs")
