require(pomp)
require(plyr)
require(reshape2)

WHO.situation.report.Oct.1 <- '
week,Guinea,Liberia,SierraLeone
1,2.244,,
2,2.244,,
3,0.073,,
4,5.717,,
5,3.954,,
6,5.444,,
7,3.274,,
8,5.762,,
9,7.615,,
10,7.615,,
11,27.392,,
12,17.387,,
13,27.115,,
14,29.29,,
15,27.84,,
16,16.345,,
17,10.917,,
18,11.959,,
19,11.959,,
20,8.657,,
21,26.537,,
22,47.764,3.517,
23,26.582,1.043,5.494
24,32.967,18,57.048
25,18.707,16.34,76.022
26,24.322,13.742,36.768
27,4.719,10.155,81.929
28,7.081,24.856,102.632
29,8.527,53.294,69.823
30,92.227,70.146,81.783
31,26.423,139.269,99.775
32,16.549,65.66,88.17
33,36.819,240.645,90.489
34,92.08,274.826,161.54
35,101.03,215.56,168.966
36,102.113,388.553,186.144
37,83.016,410.299,220.442
38,106.674,300.989,258.693
39,55.522,240.237,299.546
'

## Population sizes in Guinea, Liberia, and Sierra Leone (census 2014)
populations <- c(Guinea=10628972,Liberia=4092310,SierraLeone=6190280)
populations["WestAfrica"] <- sum(populations)

dat <- read.csv(text=WHO.situation.report.Oct.1,stringsAsFactors=FALSE)
dat <- melt(dat,id="week",variable.name="country",value.name="cases")
mutate(dat,deaths=NA) -> dat

ebolaModel <- function (country=c("Guinea", "SierraLeone", "Liberia", "WestAfrica"),
                        data = NULL,
                        timestep = 0.01, nstageE = 3L,
                        type = c("raw","cum"), na.rm = FALSE, least.sq = FALSE) {

  type <- match.arg(type)
  ctry <- match.arg(country)
  pop <- unname(populations[ctry])

  ## Incubation period is supposed to be Gamma distributed with shape parameter 3 and mean 11.4 days
  ## The discrete-time formula is used to calculate the corresponding alpha (cf He et al., Interface 2010)
  ## Case-fatality ratio is fixed at 0.7 (cf WHO Ebola response team, NEJM 2014)
  incubation_period <- 11.4/7
  infectious_period <- 7/7
  index_case <- 10/pop
  dt <- timestep

  theta <- c(N=pop,R0=1.4,
             alpha=-1/(nstageE*dt)*log(1-nstageE*dt/incubation_period),
             gamma=-log(1-dt/infectious_period)/dt,
             rho=0.2,cfr=0.7,
             k=0,
             S_0=1-index_case,E_0=index_case/2-5e-9,
             I_0=index_case/2-5e-9,R_0=1e-8)

  if (is.null(data)) {
    if (ctry=="WestAfrica") {
      dat <- ddply(dat,~week,summarize,
                   cases=sum(cases,na.rm=TRUE),
                   deaths=sum(deaths,na.rm=TRUE))
    } else {
      dat <- subset(dat,country==ctry,select=-country)
    }
  } else {
    dat <- data
  }
    
  if (na.rm) {
    dat <- mutate(subset(dat,!is.na(cases)),week=week-min(week)+1)
  }
  if (type=="cum") {
    dat <- mutate(dat,cases=cumsum(cases),deaths=cumsum(deaths))
  }

  ## Create the pomp object
  pomp(
       data=dat,
       times="week",
       t0=0,
       params=theta,
       obsnames=c("cases","deaths"),
       statenames=c("S","E1","I","R","N_EI","N_IR"),
       zeronames=if (type=="raw") c("N_EI","N_IR") else character(0),
       paramnames=c("N","R0","alpha","gamma","rho","k","cfr",
         "S_0","E_0","I_0","R_0"),
       nstageE=as.integer(nstageE),
       PACKAGE="pompExamples",
       dmeasure=if (least.sq) "_ebola_dObsLS" else "_ebola_dObs",
       rmeasure=if (least.sq) "_ebola_rObsLS" else "_ebola_rObs",
       rprocess=discrete.time.sim(step.fun="_ebola_rSim",delta.t=timestep),
       skeleton="_ebola_skel",
       skeleton.type="vectorfield",
       parameter.transform="_ebola_par_trans",
       parameter.inv.transform="_ebola_par_untrans",
       initializer=function (params, t0, nstageE, ...) {
         all.state.names <- c("S",paste0("E",1:nstageE),"I","R","N_EI","N_IR")
         comp.names <- c("S",paste0("E",1:nstageE),"I","R")
         x0 <- setNames(numeric(length(all.state.names)),all.state.names)
         frac <- c(params["S_0"],rep(params["E_0"]/nstageE,nstageE),params["I_0"],params["R_0"])
         x0[comp.names] <- round(params["N"]*frac/sum(frac))
         x0
       }
       ) -> po
}

c("ebolaModel")
