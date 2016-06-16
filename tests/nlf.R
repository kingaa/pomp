library(pomp)

set.seed(583615606L)

pompExample(ou2)
estnames=c("alpha.2","alpha.3","tau")
theta.truth <- coef(ou2)
theta.guess <- theta.truth
theta.guess[estnames] <- theta.guess[estnames]*1.5

m1 <- nlf(
    object=ou2,
    start=theta.truth,
    lags=c(4,6),
    nconverge=100,
    nasymp=2000,
    eval.only=TRUE,
    seed=426094906L,
    lql.frac = 0.025
)

m2 <- nlf(
    m1,
    est=estnames,
    maxit=500,
    method="Nelder-Mead"
)

m3 <- nlf(
    object=ou2,
    start=theta.guess,
    lags=c(4,6),
    nconverge=100,
    nasymp=2000,
    maxit=500,
    method="Nelder-Mead",
    eval.only=TRUE,
    seed=426094906L,
    lql.frac = 0.025
)

m4 <- nlf(
    m3,
    est=estnames,
    method="subplex",
    seed=300983678L
)

stopifnot(max(abs(1-c(coef(m4,estnames),se=m4$se,value=logLik(m4))/c(-0.51,0.30,1.3,0.042,0.031,0.42,-549)))<0.03)
stopifnot(max(abs(1-c(coef(m2,estnames),se=m2$se,value=logLik(m4))/c(-0.47,0.31,1.4,0.030,0.044,0.42,-550)))<0.03)

m5 <- nlf(m3,tensor=TRUE,period=10,est=estnames,seed=427458351L)
m5 <- nlf(m5,seed=469007824L)
stopifnot(max(abs(1-c(coef(m5,estnames),m5$se,logLik(m5))/c(-0.46,0.32,1.58,0.033,0.0328,0.417,-548)))<0.03)
capture.output(m6 <- nlf(m5,seed=469007824L,tensor=FALSE,verbose=TRUE)) -> msg
stopifnot(length(msg)==8)
stopifnot(sum(grepl("fitted param",msg))==6)
