library(pomp)

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
    gamma=1/3,
    Beta=1.4,Beta_sd=0,
    pop=1400,
    rho=0.9,sigma=3.6,
    S_0=0.999,I_0=0.001,R_0=0
  ),
  rprocess=euler.sim(
    step.fun=Csnippet("
      double rate[2];	  	// transition rates
      double trans[2];		// transition numbers
      double dW;
      dW = rgammawn(Beta_sd,dt); // gamma noise, mean=dt, variance=(beta_sd^2 dt)
      // compute the transition rates
      rate[0] = Beta*I/pop*dW/dt;// force of infection
      rate[1] = gamma;		       // recovery
      // compute the transition numbers
      reulermultinom(1,S,&rate[0],dt,&trans[0]);
      reulermultinom(1,I,&rate[1],dt,&trans[1]);
      // balance the equations
      S -= trans[0];
      I += trans[0]-trans[1];
      R += trans[1];
      H += trans[1];		// cases are cumulative recoveries
      if (Beta_sd > 0.0)  W += (dW-dt)/Beta_sd; // mean = 0, variance = dt"
    ),
    delta.t=1/12
  ),
  skeleton=vectorfield(Csnippet("
      double rate[2];	  	// transition rates
      double term[2];		// transition numbers
      // compute the transition rates
      rate[0] = Beta*I/pop;   // force of infection
      rate[1] = gamma;		    // recovery
      // compute the several terms
      term[0] = rate[0]*S;
      term[1] = rate[1]*I;
      // balance the equations
      DS = -term[0];
      DI = term[0]-term[1];
      DH = DR = term[1];
      DW = 0;		    	// no noise, so no noise accumulation"
  )),
  rinit=Csnippet("
    double m = pop/(S_0+I_0+R_0);
    S = nearbyint(m*S_0);
    I = nearbyint(m*I_0);
    R = nearbyint(m*R_0);
    H = 0;
    W = 0;"
  ),
  rmeasure=Csnippet("
    reports = rnbinom_mu(1/sigma/sigma,rho*H);"
  ),
  dmeasure=Csnippet("
    lik = dnbinom_mu(reports,1/sigma/sigma,rho*H,give_log);"
  ),
  statenames=c("S","I","R","H","W"),
  paramnames=c(
    "gamma","Beta","Beta_sd","pop","rho","sigma",
    "S_0","I_0","R_0"  # note we rely on these being adjacent
  ),
  zeronames=c("H"),
  partrans=parameter_trans(
    toEst=Csnippet("
      Tgamma = log(gamma);
      TBeta = log(Beta);
      TBeta_sd = log(Beta_sd);
      Trho = logit(rho);
      Tsigma = log(sigma);
      to_log_barycentric(&TS_0,&S_0,3);"
    ),
    fromEst=Csnippet("
      Tgamma = exp(gamma);
      TBeta = exp(Beta);
      TBeta_sd = exp(Beta_sd);
      Trho = expit(rho);
      Tsigma = exp(sigma);
      from_log_barycentric(&TS_0,&S_0,3);"
    )
  )
) -> bbs

c("bbs")
