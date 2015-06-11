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


paruntrans <- Csnippet('
  const double *IC = &S_0;
  double *TIC = &TS_0;
  TN = log(N);
  TR0 = log(R0);
  Talpha = log(alpha);
  Tgamma = log(gamma);
  Trho = logit(rho);
  Tk = log(k);
  Tcfr = logit(cfr);
  to_log_barycentric(TIC,IC,4);
')

partrans <- Csnippet('
  const double *IC = &S_0;
  double *TIC = &TS_0;
  TN = exp(N);
  TR0 = exp(R0);
  Talpha = exp(alpha);
  Tgamma = exp(gamma);
  Trho = expit(rho);
  Tk = exp(k);
  Tcfr = expit(cfr);
  from_log_barycentric(TIC,IC,4);
')

##  Observation model: hierarchical model for cases and deaths
## p(R_t, D_t| C_t) = p(R_t | C_t) * p(D_t | C_t, R_t)
## p(R_t | C_t): Negative binomial with mean rho * C_t and dispersion parameter 1 / k
## p(D_t | C_t, R_t): Binomial B(R_t, cfr)

dObs <- Csnippet('
  double f;
  if (k > 0.0)
    f = dnbinom_mu(nearbyint(cases),1.0/k,rho*N_EI,1);
  else
    f = dpois(nearbyint(cases),rho*N_EI,1);
  lik = (give_log) ? f : exp(f);
')

dObsLS <- Csnippet('
  double f;
  f = dnorm(cases,rho*N_EI,k,1);
  lik = (give_log) ? f : exp(f);
')

rObs <- Csnippet('
  if (k > 0) {
    cases = rnbinom_mu(1.0/k,rho*N_EI);
    deaths = rnbinom_mu(1.0/k,rho*cfr*N_IR);
  } else {
    cases = rpois(rho*N_EI);
    deaths = rpois(rho*cfr*N_IR);
  }')

rObsLS <- Csnippet('
  cases = rnorm(rho*N_EI,k);
  deaths = NA_REAL;
')

rSim <- Csnippet('
  double lambda, beta;
  double *E = &E1;
  beta = R0 * gamma; // Transmission rate
  lambda = beta * I / N; // Force of infection
  int i;

  // Transitions
  // From class S
  double transS = rbinom(S, 1.0 - exp(- lambda * dt)); // No of infections
  // From class E
  double transE[nstageE]; // No of transitions between classes E
  for(i = 0; i < nstageE; i++){
    transE[i] = rbinom(E[i], 1.0 - exp(- nstageE * alpha * dt));
  }
  // From class I
  double transI = rbinom(I, 1.0 - exp(- gamma * dt)); // No of transitions I->R

  // Balance the equations
  S -= transS;
  E[0] += transS - transE[0];
  for(i=1; i < nstageE; i++) {
    E[i] += transE[i-1] - transE[i];
  }
  I += transE[nstageE - 1] - transI;
  R += transI;
  N_EI += transE[nstageE - 1]; // No of transitions from E to I
  N_IR += transI; // No of transitions from I to R
')

skel <- Csnippet('
  double lambda, beta;
  const double *E = &E1;
  double *DE = &DE1;
  beta = R0 * gamma; // Transmission rate
  lambda = beta * I / N; // Force of infection
  int i;

  // Balance the equations
  DS = - lambda * S;
  DE[0] = lambda * S - nstageE * alpha * E[0];
  for (i=1; i < nstageE; i++)
    DE[i] = nstageE * alpha * (E[i-1]-E[i]);
  DI = nstageE * alpha * E[nstageE-1] - gamma * I;
  DR = gamma * I;
  DN_EI = nstageE * alpha * E[nstageE-1];
  DN_IR = gamma * I;
')

initlzr <- Csnippet("
  int j;
  double m, *E;
  E = &E1;
  m = N/(S_0+E_0+I_0+R_0);
  S = nearbyint(m*S_0);
  I = nearbyint(m*I_0);
  R = nearbyint(m*R_0);
  m = nearbyint(m*E_0/nstageE);
  for (j = 0; j < nstageE; j++) E[j] = m;
")


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
  nstageE <- as.integer(nstageE)

  globs <- paste0("static int nstageE = ",nstageE,";");

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
       globals=globs,
       obsnames=c("cases","deaths"),
       statenames=c("S",sprintf("E%d",1:nstageE),"I","R","N_EI","N_IR"),
       zeronames=if (type=="raw") c("N_EI","N_IR") else character(0),
       paramnames=c("N","R0","alpha","gamma","rho","k","cfr",
         "S_0","E_0","I_0","R_0"),
       dmeasure=if (least.sq) dObsLS else dObs,
       rmeasure=if (least.sq) rObsLS else rObs,
       rprocess=discrete.time.sim(step.fun=rSim,delta.t=timestep),
       skeleton=skel,
       skeleton.type="vectorfield",
       fromEstimationScale=partrans,
       toEstimationScale=paruntrans,
       initializer=initlzr
       ) -> po
}

c("ebolaModel")
