##' Ebola outbreak, West Africa
##'
##' Data and models for the 2013--2015 outbreak of Ebola virus disease in West Africa.
##'
##' The data include case counts according to the WHO situation report dated 1 October 2014.
##' The models are described in King et al. (2015).
##'
##' @name ebola
##' @docType data
##' @rdname ebola
##' @include pomp.R
##' @family datasets
##' @family pomp examples
##'
##' @references
##' A. A. King, M. Domenech de Cellès, F. M. G. Magpantay, and P. Rohani.
##' Avoidable errors in the modelling of outbreaks of emerging pathogens, with special reference to Ebola.
##' Proceedings of the Royal Society of London, Series B, 282: 20150347, 2015.
##'
##' WHO Ebola Response Team.
##' Ebola virus disease in West Africa---the first 9 months of the epidemic and forward projections.
##' New England Journal of Medicine 371: 1481--1495, 2014.
##'
##' D. He, E. L. Ionides, & A. A. King.
##' Plug-and-play inference for disease dynamics:
##' measles in large and small populations as a case study.
##' J. R. Soc. Interface, 7: 271--283, 2010.
##'
##' @section Model structure:
##' The observation model is a hierarchical model for cases and deaths:
##'   \deqn{p(R_t, D_t| C_t) = p(R_t | C_t) * p(D_t | C_t, R_t).}
##' Here, \eqn{p(R_t | C_t)} is negative binomial with mean \eqn{\rho C_t} and dispersion parameter \eqn{1/k};
##' \eqn{p(D_t | C_t, R_t)} is binomial with size \eqn{R_t} and probability equal to the case fatality rate \code{cfr}.
##'
##' The default incubation period is supposed to be Gamma distributed with shape parameter `nstageE` and mean 11.4 days.
##' The discrete-time formula is used to calculate the corresponding `alpha`` (cf. He et al. 2010).
##' The case-fatality ratio (`cfr`) is fixed at 0.7 (cf. WHO Ebola Response Team 2014).
NULL

##' @rdname ebola
##' @aliases ebolaWHO

##' @rdname ebola
##' @aliases ebolaModel
##'
##' @param country name of country or \dQuote{WestAfrica} for the regional outbreak.
##' @param data if NULL, the situation report data (see \code{ebolaWHO}) for the appropriate country or region will be used.
##' Providing a dataset here will override this behavior.
##' @param timestep duration (in weeks) of Euler time-step for the simulations.
##' @param nstageE integer; number of incubation stages.
##' @param incubation_period,infectious_period mean duration (in weeks) of the incubation and infectious periods.
##' @param type should raw or cumulative case counts be used as the data?
##' @param na.rm should missing data be removed?
##' @param least.sq logical; if \code{TRUE}, ordinary least squares will be used.
##'
##' @export
ebolaModel <- function (
  country=c("Guinea", "SierraLeone", "Liberia", "WestAfrica"),
  data = NULL,
  timestep = 0.01, nstageE = 3L,
  R0 = 1.4, rho = 0.2, cfr = 0.7, k = 0,
  incubation_period = 11.4/7, infectious_period = 7/7,
  type = c("raw","cum"), na.rm = FALSE, least.sq = FALSE
) {

  ## Population sizes in Guinea, Liberia, and Sierra Leone (census 2014)
  populations <- c(Guinea=10628972,Liberia=4092310,SierraLeone=6190280)
  populations["WestAfrica"] <- sum(populations)

  pomp2::ebolaWHO %>%
    melt(id.vars="week",variable.name="country",value.name="cases") -> dat
  dat$deaths <- NA

  type <- match.arg(type)
  ctry <- match.arg(country)
  pop <- unname(populations[ctry])
  nstageE <- as.integer(nstageE)
  dt <- as.numeric(timestep)

  index_case <- 10/pop

  globs <- paste0("static int nstageE = ",nstageE,";");

  theta <- c(
    N=pop,R0=1.4,
    rho=rho,cfr=cfr,k=k,
    alpha=-1/(nstageE*dt)*log(1-nstageE*dt/incubation_period),
    gamma=-log(1-dt/infectious_period)/dt,
    S_0=1-index_case,
    E_0=index_case/2-5e-9,
    I_0=index_case/2-5e-9,
    R_0=1e-8
  )

  if (is.null(data)) {
    if (ctry == "WestAfrica") {
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

  dObs <- Csnippet("
    double f;
    if (k > 0.0)
      f = dnbinom_mu(nearbyint(cases),1.0/k,rho*N_EI,1);
    else
      f = dpois(nearbyint(cases),rho*N_EI,1);
    lik = (give_log) ? f : exp(f);
    "
  )

  dObsLS <- Csnippet("
    double f;
    f = dnorm(cases,rho*N_EI,k,1);
    lik = (give_log) ? f : exp(f);
    "
  )

  rObs <- Csnippet("
    if (k > 0) {
      cases = rnbinom_mu(1.0/k,rho*N_EI);
      deaths = rnbinom_mu(1.0/k,rho*cfr*N_IR);
    } else {
      cases = rpois(rho*N_EI);
      deaths = rpois(rho*cfr*N_IR);
    }
    "
  )

  rObsLS <- Csnippet("
    cases = rnorm(rho*N_EI,k);
    deaths = NA_REAL;
    "
  )

  rInit <- Csnippet("
    double m = N/(S_0+E_0+I_0+R_0);
    double *E = &E1;
    int j;
    S = nearbyint(m*S_0);
    for (j = 0; j < nstageE; j++) E[j] = nearbyint(m*E_0/nstageE);
    I = nearbyint(m*I_0);
    R = nearbyint(m*R_0);
    N_EI = 0;
    N_IR = 0;
    "
  )

  rSim <- Csnippet("
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
    "
  )

  skel <- Csnippet("
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
    "
  )

  pomp(
    data=dat,
    times="week",
    t0=0,
    params=theta,
    globals=globs,
    rinit=rInit,
    rprocess=discrete_time(rSim,delta.t=timestep),
    dmeasure=if (least.sq) dObsLS else dObs,
    rmeasure=if (least.sq) rObsLS else rObs,
    skeleton=vectorfield(skel),
    partrans=parameter_trans(
      log=c("N","R0","alpha","gamma","k"),
      logit=c("rho","cfr"),
      barycentric=c("S_0","E_0","I_0","R_0")
    ),
    statenames=c("S",sprintf("E%d",1:nstageE),"I","R","N_EI","N_IR"),
    paramnames=c("N","R0","alpha","gamma","rho","k","cfr",
      "S_0","E_0","I_0","R_0"),
    accumvars=if (type=="raw") c("N_EI","N_IR") else character(0)
  )

}
