##' Model of cholera transmission for historic Bengal.
##'
##' \code{dacca} constructs a \sQuote{pomp} object containing census and cholera
##' mortality data from the Dacca district of the former British province of
##' Bengal over the years 1891 to 1940 together with a stochastic differential
##' equation transmission model.
##' The model is that of King et al. (2008).
##' The parameters are the MLE for the SIRS model with seasonal reservoir.
##'
##' Data are provided courtesy of Dr. Menno J. Bouma, London School of Tropical
##' Medicine and Hygiene.
##'
##' @return
##' \code{dacca} returns a \sQuote{pomp} object containing the model, data, and MLE
##' parameters, as estimated by King et al. (2008).
##'
##' @name dacca
##' @docType data
##' @family pomp examples
##' @family pomp datasets
##' @include pomp.R
##' @importFrom stats smooth.spline predict
##'
##' @references
##'
##' \King2008
##'
##' @keywords models
##' @example examples/dacca.R
##'
NULL

## pomp object encoding the "SIRS model with seasonal reservoir" of
##   King, A. A., Ionides, E. L., Pascual, M., & Bouma, M. J.
##   Inapparent infections and cholera dynamics.
##   Nature 454:877-880 (2008)
## Data are cholera deaths and decadal census figures from the Dacca district of Bengal province, 1891-1941.

##' @rdname dacca
##'
##' @param gamma recovery rate
##' @param eps rate of waning of immunity for severe infections
##' @param rho rate of waning of immunity for inapparent infections
##' @param delta baseline mortality rate
##' @param deltaI cholera mortality rate
##' @param clin fraction of infections that lead to severe infection
##' @param alpha transmission function exponent
##' @param beta_trend slope of secular trend in transmission
##' @param logbeta seasonal transmission rates
##' @param sd_beta environmental noise intensity
##' @param tau measurement error s.d.
##' @param logomega seasonal environmental reservoir parameters
##' @param S_0 initial susceptible fraction
##' @param I_0 initial fraction of population infected
##' @param Y_0 initial fraction of the population in the Y class
##' @param R1_0,R2_0,R3_0 initial fractions in the respective R classes
##'
##' @export
dacca <- function (
  gamma = 20.8, eps = 19.1, rho = 0,
  delta = 0.02, deltaI = 0.06, clin = 1, alpha = 1,
  beta_trend = -0.00498,
  logbeta = c(0.747, 6.38, -3.44, 4.23, 3.33, 4.55),
  logomega = log(c(0.184, 0.0786, 0.0584, 0.00917, 0.000208, 0.0124)),
  sd_beta = 3.13, tau = 0.23,
  S_0 = 0.621, I_0 = 0.378, Y_0 = 0,
  R1_0 = 0.000843, R2_0 = 0.000972, R3_0 = 1.16e-07
)
{
  nrstage <- 3L
  nbasis <- length(logbeta)
  if (length(logomega) != nbasis)
    pStop(sQuote("logbeta")," and ",sQuote("logomega")," should be of equal length.")

  mle <- c(
    gamma=gamma,eps=eps,rho=rho,
    delta=delta, deltaI=deltaI, clin=clin, alpha=alpha,
    beta_trend=beta_trend,
    logbeta=unname(logbeta),
    logomega=unname(logomega),
    sd_beta=sd_beta, tau=tau,
    S_0=S_0, I_0=I_0, Y_0=Y_0, R1_0=R1_0, R2_0=R2_0, R3_0=R3_0
  )

  census <- data.frame(
    year = c(1891L, 1901L, 1911L, 1921L, 1931L, 1941L),
    census = c(2420656L, 2649522L, 2960402L, 3125967L, 3432577L, 4222142L)
  )

  cholera <- data.frame(
    time=seq(from=1891+1/12,to=1941,by=1/12),
    cholera.deaths = c(
      2641L, 939L, 905L, 1219L, 368L, 78L, 29L, 12L, 30L, 44L, 270L, 1149L,
      633L, 501L, 855L, 1271L, 666L, 101L, 62L, 23L, 20L, 28L, 461L,
      892L, 751L, 170L, 253L, 906L, 700L, 98L, 57L, 72L, 471L, 4217L,
      5168L, 4747L, 2380L, 852L, 1166L, 2122L, 576L, 60L, 53L, 62L,
      241L, 403L, 551L, 739L, 862L, 348L, 490L, 5596L, 1180L, 142L,
      41L, 28L, 39L, 748L, 3934L, 3562L, 587L, 311L, 1639L, 1903L,
      601L, 110L, 32L, 19L, 82L, 420L, 1014L, 1073L, 416L, 168L, 909L,
      1355L, 447L, 59L, 13L, 21L, 43L, 109L, 338L, 470L, 489L, 394L,
      483L, 842L, 356L, 29L, 17L, 16L, 57L, 110L, 488L, 1727L, 1253L,
      359L, 245L, 549L, 215L, 9L, 7L, 31L, 236L, 279L, 819L, 1728L,
      1942L, 1251L, 3521L, 3412L, 290L, 46L, 35L, 14L, 79L, 852L, 2951L,
      2656L, 607L, 172L, 325L, 2191L, 584L, 58L, 38L, 8L, 22L, 50L,
      380L, 2059L, 938L, 389L, 767L, 1882L, 286L, 94L, 61L, 10L, 106L,
      281L, 357L, 1388L, 810L, 306L, 381L, 1308L, 702L, 87L, 9L, 14L,
      36L, 46L, 553L, 1302L, 618L, 147L, 414L, 768L, 373L, 39L, 10L,
      36L, 151L, 1130L, 3437L, 4041L, 1415L, 207L, 92L, 128L, 147L,
      32L, 7L, 59L, 426L, 2644L, 2891L, 4249L, 2291L, 797L, 680L, 1036L,
      404L, 41L, 19L, 12L, 10L, 121L, 931L, 2158L, 1886L, 803L, 397L,
      613L, 132L, 48L, 17L, 22L, 26L, 34L, 344L, 657L, 117L, 75L, 443L,
      972L, 646L, 107L, 18L, 6L, 9L, 5L, 12L, 142L, 133L, 189L, 1715L,
      3115L, 1412L, 182L, 50L, 37L, 77L, 475L, 1730L, 1489L, 620L,
      190L, 571L, 1558L, 440L, 27L, 7L, 14L, 93L, 1462L, 2467L, 1703L,
      1262L, 458L, 453L, 717L, 232L, 26L, 16L, 18L, 9L, 78L, 353L,
      897L, 777L, 404L, 799L, 2067L, 613L, 98L, 19L, 26L, 47L, 171L,
      767L, 1896L, 887L, 325L, 816L, 1653L, 355L, 85L, 54L, 88L, 609L,
      882L, 1363L, 2178L, 580L, 396L, 1493L, 2154L, 683L, 78L, 19L,
      10L, 27L, 88L, 1178L, 1862L, 611L, 478L, 2697L, 3395L, 520L,
      67L, 41L, 36L, 209L, 559L, 971L, 2144L, 1099L, 494L, 586L, 508L,
      269L, 27L, 19L, 21L, 12L, 22L, 333L, 676L, 487L, 262L, 535L,
      979L, 170L, 25L, 9L, 19L, 13L, 45L, 229L, 673L, 432L, 107L, 373L,
      1126L, 339L, 19L, 11L, 3L, 15L, 101L, 539L, 709L, 200L, 208L,
      926L, 1783L, 831L, 103L, 37L, 17L, 33L, 179L, 426L, 795L, 481L,
      491L, 773L, 936L, 325L, 101L, 22L, 25L, 24L, 88L, 633L, 513L,
      298L, 93L, 687L, 1750L, 356L, 33L, 2L, 18L, 70L, 648L, 2471L,
      1270L, 616L, 193L, 706L, 1372L, 668L, 107L, 58L, 21L, 23L, 93L,
      318L, 867L, 332L, 118L, 437L, 2233L, 491L, 27L, 7L, 21L, 96L,
      360L, 783L, 1492L, 550L, 176L, 633L, 922L, 267L, 91L, 42L, 4L,
      10L, 7L, 43L, 377L, 563L, 284L, 298L, 625L, 131L, 35L, 12L, 8L,
      9L, 83L, 502L, 551L, 256L, 198L, 664L, 1701L, 425L, 76L, 17L,
      9L, 16L, 5L, 141L, 806L, 1603L, 587L, 530L, 771L, 511L, 97L,
      35L, 39L, 156L, 1097L, 1233L, 1418L, 1125L, 420L, 1592L, 4169L,
      1535L, 371L, 139L, 55L, 85L, 538L, 1676L, 1435L, 804L, 370L,
      477L, 394L, 306L, 132L, 84L, 87L, 53L, 391L, 1541L, 1859L, 894L,
      326L, 853L, 1891L, 1009L, 131L, 77L, 63L, 66L, 33L, 178L, 1003L,
      1051L, 488L, 911L, 1806L, 837L, 280L, 132L, 76L, 381L, 1328L,
      2639L, 2164L, 1082L, 326L, 254L, 258L, 119L, 106L, 93L, 29L,
      17L, 17L, 17L, 46L, 79L, 135L, 1290L, 2240L, 561L, 116L, 24L,
      15L, 33L, 18L, 16L, 38L, 26L, 45L, 151L, 168L, 57L, 32L, 29L,
      27L, 20L, 106L, 1522L, 2013L, 434L, 205L, 528L, 634L, 195L, 45L,
      33L, 19L, 20L, 46L, 107L, 725L, 572L, 183L, 2199L, 4018L, 428L,
      67L, 31L, 8L, 44L, 484L, 1324L, 2054L, 467L, 216L, 673L, 887L,
      353L, 73L, 46L, 15L, 20L, 27L, 25L, 38L, 158L, 312L, 1226L, 1021L,
      222L, 90L, 31L, 93L, 368L, 657L, 2208L, 2178L, 702L, 157L, 317L,
      146L, 63L, 27L, 22L, 23L, 28L, 225L, 483L, 319L, 120L, 59L, 274L,
      282L, 155L, 31L, 16L, 15L, 12L, 14L, 14L, 42L
    )
  )

  rinit <- Csnippet("
    int k;
    double sum = S_0+I_0+Y_0;
    double *R = &R1;
    const double *R0 = &R1_0;
    for (k = 0; k < nrstage; k++) sum += R0[k];
    S = nearbyint(pop*S_0/sum);
    I = nearbyint(pop*I_0/sum);
    Y = nearbyint(pop*Y_0/sum);
    for (k = 0; k < nrstage; k++) R[k] = nearbyint(pop*R0[k]/sum);
    W = 0;
    deaths = 0;
    count = 0;
  ")

  norm_rmeas <- Csnippet("
    double v, tol = 1.0e-18;
    v = deaths*tau;
    if ((count > 0) || (!(R_FINITE(v)))) {
    cholera_deaths = R_NaReal;
    } else {
    cholera_deaths = rnorm(deaths,v+tol);
    }
  ")

  norm_dmeas <- Csnippet("
    double v, tol = 1.0e-18;
    v = deaths*tau;
    if ((count>0.0) || (!(R_FINITE(v)))) {
    lik = tol;
    } else {
    lik = dnorm(cholera_deaths,deaths,v+tol,0)+tol;
    }
    if (give_log) lik = log(lik);
  ")

  ## two-path SIRS cholera model using SDEs
  ## exponent (alpha) on I/n
  ## only "severe" infections are infectious
  ## truncation is not used
  ## instead, particles with negative states are killed

  cholmodel_one <- Csnippet("
    double births;
    double infections;
    double sdeaths;
    double ideaths;
    double ydeaths;
    double rdeaths[nrstage];
    double disease;
    double wanings;
    double passages[nrstage+1];
    double effI;
    double neps;
    double beta;
    double omega;
    double dw;
    double *pt;
    int j;

    if (count != 0.0) return;

    neps = eps*nrstage;

    beta = exp(dot_product(nbasis,&seas_1,&logbeta1)+beta_trend*trend);
    omega = exp(dot_product(nbasis,&seas_1,&logomega1));

    dw = rnorm(0,sqrt(dt));     // white noise

    effI = pow(I/pop,alpha);
    births = dpopdt + delta*pop;        // births

    passages[0] = gamma*I;      // recovery
    ideaths = delta*I;          // natural i deaths
    disease = deltaI*I;         // disease death
    ydeaths = delta*Y;          // natural rs deaths
    wanings = rho*Y;            // loss of immunity

    for (pt = &R1, j = 0; j < nrstage; j++, pt++) {
    rdeaths[j] = *pt*delta;     // natural R deaths
    passages[j+1] = *pt*neps;   // passage to the next immunity class
    }

    infections = (omega+(beta+sd_beta*dw/dt)*effI)*S; // infection
    sdeaths = delta*S;          // natural S deaths

    S += (births - infections - sdeaths + passages[nrstage] + wanings)*dt;
    I += (clin*infections - disease - ideaths - passages[0])*dt;
    Y += ((1-clin)*infections - ydeaths - wanings)*dt;
    for (pt = &R1, j = 0; j < nrstage; j++, pt++)
    *pt += (passages[j] - passages[j+1] - rdeaths[j])*dt;
    deaths += disease*dt;               // cumulative deaths due to disease
    W += dw;

    // check for violations of positivity constraints
    // nonzero 'count' variable signals violation
    if (S < 0.0) {
    S = 0.0; I = 0.0; Y = 0.0;
    count += 1;
    }
    if (I < 0.0) {
    I = 0.0; S = 0.0;
    count += 1e3;
    }
    if (Y < 0.0) {
    Y = 0.0; S = 0.0;
    count += 1e6;
    }
    if (deaths < 0.0) {
    deaths = 0.0;
    count += 1e9;
    }
    for (pt = &R1, j = 0; j < nrstage-1; j++, pt++) {
    if (*pt < 0.0) {
    *pt = 0.0; *(pt+1) = 0.0;
    count += 1e12;
    }
    }
    if (*pt < 0.0) {
    *pt = 0.0; S = 0.0;
    count += 1e12;
    }
    "
  )

  t0 <- with(cholera,2*time[1]-time[2])

  pomp(
    data=cholera,
    times='time',
    t0=t0,
    params=mle,
    cfile="dacca_source",
    globals = sprintf("int nrstage = %d, nbasis = %d;",nrstage,nbasis),
    rprocess = euler(
      step.fun = cholmodel_one,
      delta.t=1/240
    ),
    dmeasure = norm_dmeas,
    rmeasure=norm_rmeas,
    partrans=parameter_trans(
      log=c("tau","gamma","eps","delta","deltaI","sd_beta","alpha","rho"),
      logit="clin",
      barycentric=c("S_0","I_0","Y_0",sprintf("R%01d_0",1:nrstage))
    ),
    rinit=rinit,
    covar=covariate_table(
      t=seq(from=t0,to=max(cholera$time)+2/12,by=0.01),
      seas=periodic_bspline_basis(t-1/12,nbasis=nbasis,degree=3,period=1),
      pop=predict(smooth.spline(x=census$year,y=census$census),x=t)$y,
      dpopdt=predict(smooth.spline(x=census$year,y=census$census),x=t,deriv=1)$y,
      trend=t-mean(t),
      times="t"
    ),
    accumvars = c("deaths","count"),
    statenames = c("S","I","Y",sprintf("R%d",seq_len(nrstage)),"deaths","W","count"),
    paramnames = c("tau","gamma","eps","delta","deltaI",
      "logomega1","sd_beta","beta_trend","logbeta1",
      "alpha","rho","clin","S_0","I_0","Y_0","R1_0")
  )
}
