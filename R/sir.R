##' Compartmental epidemiological models
##'
##' Simple SIR-type models implemented in various ways.
##'
##' \code{sir()} producees a \sQuote{pomp} object encoding a simple seasonal SIR model with simulated data.
##' Simulation is performed using an Euler multinomial approximation.
##'
##' \code{sir2()} has the same model implemented using Gillespie's algorithm.
##'
##' This and similar examples are discussed and constructed in tutorials
##' available on the \href{https://kingaa.github.io/pomp/}{package website}.
##'
##' @name sir_models
##' @rdname sir
##' @aliases sir sir2
##' @docType data
##' @keywords datasets models
##' @family pomp examples
##'
##' @return
##' These functions return \sQuote{pomp} objects containing simulated data.
##'
##' @example examples/sir.R
##'
NULL

##' @rdname sir
##'
##' @param gamma recovery rate
##' @param mu death rate (assumed equal to the birth rate)
##' @param iota infection import rate
##' @param rho reporting efficiency
##' @param pop overall host population size
##' @param S_0,I_0,R_0 the fractions of the host population that are susceptible, infectious, and recovered, respectively, at time zero.
##' @param beta1,beta2,beta3 seasonal contact rates
##' @param beta_sd environmental noise intensity
##' @param t0 zero time
##' @param times observation times
##' @param seed seed of the random number generator
##' @param delta.t Euler step size
##'
##' @export
sir <- function (
  gamma = 26, mu = 0.02, iota = 0.01,
  beta1 = 400, beta2 = 480, beta3 = 320,
  beta_sd = 1e-3, rho = 0.6,
  pop = 2.1e6,
  S_0 = 26/400, I_0 = 0.001, R_0 = 1-S_0-I_0,
  t0 = 0,
  times = seq(from = t0 + 1/52, to = t0 + 4, by = 1/52),
  seed=329343545, delta.t = 1/52/20
) {
  
  tt0 <- t0
  tt <- times

  simulate(
    times=tt,t0=tt0,
    seed=seed,
    params=c(
      gamma=gamma,mu=mu,iota=iota,
      beta1=beta1,beta2=beta2,beta3=beta3,
      beta_sd=beta_sd,rho=rho,
      pop=pop,
      S_0=S_0,I_0=I_0,R_0=R_0
    ),
    covar=covariate_table(
      t=seq(from=tt0,to=max(tt)+0.2,by=0.01),
      seas=periodic.bspline.basis(
        x=t,
        period=1,
        nbasis=3,
        degree=3
      ),
      times="t"
    ),
    cfile="sir_source",
    globals=Csnippet("
      static int nbasis = 3;"
    ),
    rinit=Csnippet("
      double m = pop/(S_0+I_0+R_0);
      S = nearbyint(m*S_0);
      I = nearbyint(m*I_0);
      R = nearbyint(m*R_0);
      cases = 0;
      W = 0;"
    ),
    rprocess=euler(
      step.fun=Csnippet("
        int nrate = 6;
        double rate[nrate];		  // transition rates
        double trans[nrate];		// transition numbers
        double beta;
        double dW;

        beta = dot_product(nbasis,&beta1,&seas_1);

        // gamma noise, mean=dt, variance=(beta_sd^2 dt)
        dW = rgammawn(beta_sd,dt);

        // compute the transition rates
        rate[0] = mu*pop;		// birth into susceptible class
        rate[1] = (iota+beta*I*dW/dt)/pop; // force of infection
        rate[2] = mu;			// death from susceptible class
        rate[3] = gamma;	// recovery
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
        cases += trans[3];		// cases are cumulative recoveries
        if (beta_sd > 0.0)  W += (dW-dt)/beta_sd;"
      ),
      delta.t=delta.t
    ),
    skeleton=vectorfield(Csnippet("
      int nrate = 6;
      double rate[nrate];		  // transition rates
      double term[nrate];		// terms in the equations
      double beta;
      double dW;

      beta = dot_product(nbasis,&beta1,&seas_1);

      // compute the transition rates
      rate[0] = mu*pop;		// birth into susceptible class
      rate[1] = (iota+beta*I)/pop; // force of infection
      rate[2] = mu;			// death from susceptible class
      rate[3] = gamma;	// recovery
      rate[4] = mu;			// death from infectious class
      rate[5] = mu; 		// death from recovered class

      // compute the several terms
      term[0] = rate[0];
      term[1] = rate[1]*S;
      term[2] = rate[2]*S;
      term[3] = rate[3]*I;
      term[4] = rate[4]*I;
      term[5] = rate[5]*R;

      // balance the equations
      DS = term[0]-term[1]-term[2];
      DI = term[1]-term[3]-term[4];
      DR = term[3]-term[5];
      Dcases = term[3];		// accumulate the new I->R transitions
      DW = 0;			// no noise, so no noise accumulation"
      )
    ),
    rmeasure=Csnippet("
      double mean, sd;
      double rep;
      mean = cases*rho;
      sd = sqrt(cases*rho*(1-rho));
      rep = nearbyint(rnorm(mean,sd));
      reports = (rep > 0) ? rep : 0;"
    ),
    dmeasure=Csnippet("
      double mean, sd;
      double f;
      mean = cases*rho;
      sd = sqrt(cases*rho*(1-rho));
      if (reports > 0) {
      f = pnorm(reports+0.5,mean,sd,1,0)-pnorm(reports-0.5,mean,sd,1,0);
      } else {
      f = pnorm(reports+0.5,mean,sd,1,0);
      }
      lik = (give_log) ? log(f) : f;"
    ),
    partrans=parameter_trans(
      fromEst=Csnippet("
        int k;
        const double *TBETA = &T_beta1;
        double *BETA = &beta1;
        for (k = 0; k < nbasis; k++) BETA[k] = exp(TBETA[k]);"
      ),
      toEst=Csnippet("
        int k;
        const double *BETA = &beta1;
        double *TBETA = &T_beta1;
        for (k = 0; k < nbasis; k++) TBETA[k] = log(BETA[k]);"
      ),
      log=c("gamma","mu","iota","beta_sd"),
      logit="rho",
      barycentric=c("S_0","I_0","R_0")
    ),
    statenames=c("S","I","R","cases","W"),
    obsnames="reports",
    paramnames=c(
      "gamma","mu","iota",
      "beta1","beta_sd","pop","rho",
      "S_0","I_0","R_0"
    ),
    accumvars=c("cases")
  )
}

##' @name sir2
##' @rdname sir
##' @export
sir2 <- function (
  gamma = 24, mu = 1/70, iota = 0.1,
  beta1 = 330, beta2 = 410, beta3 = 490,
  rho = 0.1, pop = 1e6,
  S_0 = 0.05, I_0 = 0.0001, R_0 = 1 - S_0 - I_0,
  t0 = 0,
  times = seq(from = t0 + 1/12, to = t0 + 10, by=1/12),
  seed=1772464524
) {

  tt0 <- t0
  tt <- times

  simulate(
    times=tt,t0=tt0,
    seed=seed,
    params=c(
      gamma=gamma,mu=mu,iota=iota,rho=rho,
      beta1=beta1,beta2=beta2,beta3=beta3,
      S_0=S_0,I_0=I_0,R_0=R_0,
      pop=pop
    ),
    cfile="sir2_source",
    covar=covariate_table(
      t=seq(tt0,max(tt)+0.2,by=0.01),
      seas=periodic.bspline.basis(
        x=t,
        period=1,
        nbasis=3,
        degree=3
      ),
      times="t"
    ),
    globals=Csnippet("static int nbasis = 3;"),
    rprocess=gillespie_hl(
      .pre="double beta;",
      birth=list(
        "rate = mu*pop;",
        c(S=1,I=0,R=0,N=1,cases=0)),
      susc.death=list(
        "rate = mu*S;",
        c(S=-1,I=0,R=0,N=-1,cases=0)),
      infection=list("
      beta = dot_product(nbasis,&beta1,&seas_1);
      rate = (beta*I+iota)*S/pop;",
        c(S=-1,I=1,N=0,R=0,cases=0)),
      inf.death=list(
        "rate = mu*I;",
        c(S=0,I=-1,R=0,N=-1,cases=0)),
      recovery=list(
        "rate = gamma*I;",
        c(S=0,I=-1,R=1,N=0,cases=1)),
      recov.death=list(
        "rate = mu*R;",
        c(S=0,I=0,R=-1,N=-1,cases=0)),
      hmax=0.05),
    skeleton=vectorfield(
      Csnippet("
        int nrate = 6;
        double rate[nrate];
        double term[nrate];
        double beta;

        beta = dot_product(nbasis,&beta1,&seas_1);

        rate[0] = mu*pop;
        rate[1] = (iota+beta*I)/pop;
        rate[2] = mu;
        rate[3] = gamma;
        rate[4] = mu;
        rate[5] = mu;

        term[0] = rate[0];
        term[1] = rate[1]*S;
        term[2] = rate[2]*S;
        term[3] = rate[3]*I;
        term[4] = rate[4]*I;
        term[5] = rate[5]*R;

        DS = term[0]-term[1]-term[2];
        DI = term[1]-term[3]-term[4];
        DR = term[3]-term[5];
        DN = term[0]-term[2]-term[4]-term[5];
        Dcases = term[3];"
      )
    ),
    rmeasure=Csnippet("
      reports = rbinom(cases,rho);"
    ),
    dmeasure=Csnippet("
      lik = dbinom(reports,cases,rho,give_log);"
    ),
    statenames=c("S","I","R","N","cases"),
    obsnames="reports",
    paramnames=c(
      "gamma","mu","iota","pop","rho",
      "beta1","beta2","beta3",
      "S_0","I_0","R_0"
    ),
    accumvars=c("cases"),
    partrans=parameter_trans(
      log=c("gamma","mu","iota",sprintf("beta%d",1:3)),
      logit="rho",
      barycentric=c("S_0","I_0","R_0")
    ),
    rinit=Csnippet("
      double m;
      m = pop/(S_0+I_0+R_0);
      S = nearbyint(m*S_0);
      I = nearbyint(m*I_0);
      N = nearbyint(pop);
      R = nearbyint(m*R_0);
      cases = 0;"
    )
  )

}
