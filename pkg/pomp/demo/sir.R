require(pomp)

## negative binomial measurement model
dmeas <- "
  double prob = theta/(theta+rho*incid);
  lik = dnbinom(cases,theta,prob,give_log);
"
rmeas <- "
  double prob = theta/(theta+rho*incid);
  cases = rnbinom(theta,prob);
"
## SIR process model with extra-demographic stochasticity
step.fn <- '
  int nrate = 6;
  double rate[nrate];		// transition rates
  double trans[nrate];		// transition numbers
  double beta;			// transmission rate
  double dW;			// white noise increment
  int k;

  beta = beta1*seas1+beta2*seas2+beta3*seas3;

  // compute the environmental stochasticity
  dW = rgammawn(beta_sd,dt);

  // compute the transition rates
  rate[0] = mu*popsize;		// birth into susceptible class
  rate[1] = (iota+beta*I*dW/dt)/popsize; // force of infection
  rate[2] = mu;			// death from susceptible class
  rate[3] = gamma;		// recovery
  rate[4] = mu;			// death from infectious class
  rate[5] = mu; 		// death from recovered class

  // compute the transition numbers
  trans[0] = rpois(rate[0]*dt);	// births are Poisson
  reulermultinom(2,S,&rate[1],dt,&trans[1]);
  reulermultinom(2,I,&rate[3],dt,&trans[3]);
  reulermultinom(1,R,&rate[5],dt,&trans[5]);

  // balance the equations
  S += trans[0]-trans[1]-trans[2];
  I += trans[1]-trans[3]-trans[4];
  R += trans[3]-trans[5];
  incid += trans[3];		// cases are cumulative recoveries
  if (beta_sd > 0.0) W += (dW-dt)/beta_sd; // increment has mean = 0, variance = dt
'
skel <- '
  int nrate = 6;
  double rate[nrate];		// transition rates
  double term[nrate];		// transition numbers
  double beta;			// transmission rate
  double dW;			// white noise increment
  int k;
  
  beta = beta1*seas1+beta2*seas2+beta3*seas3;

  // compute the transition rates
  rate[0] = mu*popsize;		// birth into susceptible class
  rate[1] = (iota+beta*I)/popsize; // force of infection
  rate[2] = mu;			// death from susceptible class
  rate[3] = gamma;		// recovery
  rate[4] = mu;			// death from infectious class
  rate[5] = mu; 		// death from recovered class

  // compute the several terms
  term[0] = rate[0];
  term[1] = rate[1]*S;
  term[2] = rate[2]*S;
  term[3] = rate[3]*I;
  term[4] = rate[4]*I;
  term[5] = rate[5]*R;

  // assemble the differential equations
  DS = term[0]-term[1]-term[2];
  DI = term[1]-term[3]-term[4];
  DR = term[3]-term[5];
  Dincid = term[3];		// accumulate the new I->R transitions
  DW = 0;
'
## parameter transformations
partrans <- "
  double sum;
  Tgamma = exp(gamma);
  Tmu = exp(mu);
  Tiota = exp(iota);
  Tbeta1 = exp(beta1);
  Tbeta2 = exp(beta2);
  Tbeta3 = exp(beta3);
  Tbeta_sd = exp(beta_sd);
  Trho = expit(rho);
  Ttheta = exp(theta);
  TS_0 = exp(S_0);
  TI_0 = exp(I_0);
  TR_0 = exp(R_0);
  sum = TS_0+TI_0+TR_0;
  TS_0 /= sum;
  TI_0 /= sum;
  TR_0 /= sum;
"
paruntrans <- "
  double sum;
  Tgamma = log(gamma);
  Tmu = log(mu);
  Tiota = log(iota);
  Tbeta1 = log(beta1);
  Tbeta2 = log(beta2);
  Tbeta3 = log(beta3);
  Tbeta_sd = log(beta_sd);
  Trho = logit(rho);
  Ttheta = log(theta);
  sum = S_0+I_0+R_0;
  TS_0 = log(S_0/sum);
  TI_0 = log(I_0/sum);
  TR_0 = log(R_0/sum);
"

data(LondonYorke)

cbind(
      time=seq(from=1928,to=1934,by=0.01),
      as.data.frame(
                    periodic.bspline.basis(
                                           x=seq(from=1928,to=1934,by=0.01),
                                           nbasis=3,
                                           degree=3,
                                           period=1,
                                           names="seas%d"
                                           )
                    )
      ) -> covar

pompBuilder(
            name="SIR",
            data=subset(
              LondonYorke,
              subset=town=="New York"&disease=="measles"&year>=1928&year<=1933,
              select=c(time,cases)
            ),
            times="time",
            t0=1928,
            dmeasure=dmeas,
            rmeasure=rmeas,
            step.fn=step.fn,
            step.fn.delta.t=1/52/20,
            skeleton.type="vectorfield",
            skeleton=skel,
            covar=covar,
            tcovar="time",
            parameter.transform=partrans,
            parameter.inv.transform=paruntrans,
            statenames=c("S","I","R","incid","W"),
            paramnames=c(
              "gamma","mu","iota","beta1","beta2","beta3","beta.sd",
              "popsize","rho","theta","S.0","I.0","R.0"
              ), 
            zeronames=c("incid","W"),
            comp.names=c("S","I","R"),
            ic.names=c("S.0","I.0","R.0"),
            initializer=function(params, t0, comp.names, ic.names, ...) {
              x0 <- numeric(5)
              names(x0) <- c("S","I","R","incid","W")
              fracs <- params[ic.names]
              x0[comp.names] <- round(params['popsize']*fracs/sum(fracs))
              x0
            }
            ) -> po

coef(po) <- c(
              gamma=26,mu=0.02,iota=0.01,
              beta1=120,beta2=140,beta3=100,
              beta.sd=0.01,
              popsize=5e6,
              rho=0.2,theta=0.001,
              S.0=0.22,I.0=0.0018,R.0=0.78
              )

## compute a trajectory of the deterministic skeleton
tic <- Sys.time()
X <- trajectory(po,hmax=1/52,as.data.frame=TRUE)
toc <- Sys.time()
print(toc-tic)

plot(incid~time,data=X,type='l')

## simulate from the model
tic <- Sys.time()
x <- simulate(po,nsim=3,as.data.frame=TRUE)
toc <- Sys.time()
print(toc-tic)

plot(incid~time,data=x,col=as.factor(x$sim),pch=16)

coef(po) <- coef(
                 traj.match(
                            pomp(
                                 window(po,end=1930)
                                 ## window(po,end=1930),
                                 ## measurement.model=cases~norm(mean=rho*incid,sd=100)
                                 ),
                            est=c("S.0","I.0","R.0"),
                            transform=TRUE
                            )
                 )

pf <- pfilter(po,Np=1000,max.fail=100)
print(round(logLik(pf),1))

dyn.unload("SIR.so")
