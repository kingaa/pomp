capture.output({

    library(pomp)

    options(digits=3,verbose=TRUE,keep.source=TRUE)

    set.seed(398585L)
    pompExample(ou2)

    time(ou2) <- 1:20

    Np <- 100

    prior.bounds <- rbind(
        alpha.2=c(-0.55,-0.45),
        alpha.3=c(0.25,0.35)
    )

    pomp(ou2,
         rprior=function(params,...) {
             params[c("alpha.2","alpha.3")] <-
                 runif(n=1,min=prior.bounds[,1],max=prior.bounds[,2])
             params
         },
         dprior=function(params,log,...) {
             f <- sum(dunif(params,min=coef(ou2)-1,max=coef(ou2)+1,log=TRUE))
             if (log) f else exp(f)
         }
         ) -> ou2

    f1 <- bsmc(ou2,est="alpha.2",Np=100,smooth=0.02)
    try(bsmc(ou2,est="alpha.2",Np=2,smooth=0.02))
    f1 <- bsmc2(ou2,est="alpha.2",Np=100,smooth=0.02)
    try(bsmc(ou2,est=c("alpha.1","alpha.2"),Np=1,smooth=0.02))
    try(bsmc2(ou2,est="alpha.2",Np=1,smooth=0.02))
    try(bsmc2(ou2,est="alpha.2",Np=2,smooth=0.02))
    f1 <- ou2
    f1@data[,c(3,20)] <- c(10000,10000)
    try(f1 <- bsmc2(f1,est=c("alpha.2","alpha.4"),Np=100,smooth=0.01,max.fail=3))
    prop <- mvn.diag.rw(c(alpha.2=0.001,alpha.3=0.001))
    f2 <- pmcmc(ou2,Nmcmc=20,proposal=prop,Np=100)
    f3 <- ou2
    f3@data[,20] <- c(1000,1000)
    timezero(f3) <- 1
    f3 <- pfilter(f3,Np=10,filter.traj=TRUE)
    f3 <- pfilter(ou2,Np=100)
    f4 <- mif(f3,Nmif=10,rw.sd=c(alpha.2=0.01,alpha.3=0.01),cooling.fraction.50=0.1)
    f5 <- mif2(f3,Nmif=10,rw.sd=rw.sd(alpha.2=0.01,alpha.3=0.01),
               cooling.fraction.50=0.1)
    plist <- list(
        y1.mean=probe.mean(var="y1"),
        y2.mean=probe.mean(var="y2"),
        probe.acf(var="y1",lags=c(0,5)),
        probe.acf(var="y2",lags=c(0,5)),
        probe.ccf(vars=c("y1","y2"),lags=0)
    )
    f6 <- probe(ou2,probes=plist,nsim=200)
    f7 <- probe.match(f6,est=c("alpha.2","alpha.3"))
    f8 <- abc(f7,Nabc=20,est=c("alpha.2","alpha.3"),
              proposal=prop,scale=1,epsilon=20)
    f9 <- nlf(ou2,lags=c(1,2),est=c("alpha.2","alpha.3","tau"),
              nconverge=100,nasymp=2000,lql.frac=0.025,
              seed=426094906L,maxit=200,method="Nelder-Mead")
    f10 <- traj.match(f9,est=c("alpha.2","alpha.3","tau"))

    pompExample(ricker)
    try(pomp(ricker,rmeasure=Csnippet("y=rpois(N)"),statenames="N") -> po)
}) -> out
length(out)
stopifnot(sum(grepl("mif2 pfilter",out))==40)
stopifnot(sum(grepl("model codes written",out))==1)
stopifnot(sum(grepl("fitted param",out))==6)
stopifnot(sum(grepl("ABC iteration",out))==5)
stopifnot(sum(grepl("acceptance ratio:",out))==24)
stopifnot(sum(grepl("pfilter timestep",out))==88)
stopifnot(sum(grepl("mif iteration",out))==10)
stopifnot(sum(grepl("prior.mean",out))==78)
stopifnot(sum(grepl("effective sample size",out))==74)
stopifnot(sum(grepl("mif2 iteration",out))==10)
stopifnot(sum(grepl("h in",out))==1)

invisible(capture.output(pompExample(ricker)))
capture.output(simulate(pomp(ricker,rmeasure=Csnippet("y=rpois(N);"),statenames="N",
  cfile="bob",verbose=TRUE),verbose=TRUE) -> po) -> out
gsub("(\\w+)\\s.+","\\1",out,perl=TRUE)

capture.output(invisible(mif2(window(ricker,end=10),Nmif=1,Np=1,rw.sd=rw.sd(r=1),
  transform=TRUE,cooling.fraction.50=1,verbose=TRUE)),
  type="message") -> out
stopifnot(sum(grepl("filtering failure at time",out))==5)
