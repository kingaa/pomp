\donttest{
  ## A simple SIR model.

  ewmeas |>
    subset(time < 1952) |>
    pomp(
      times="time",t0=1948,
      rprocess=euler(
        Csnippet("
        int nrate = 6;
        double rate[nrate];     // transition rates
        double trans[nrate];    // transition numbers
        double dW;

        // gamma noise, mean=dt, variance=(sigma^2 dt)
        dW = rgammawn(sigma,dt);

        // compute the transition rates
        rate[0] = mu*pop;       // birth into susceptible class
        rate[1] = (iota+Beta*I*dW/dt)/pop; // force of infection
        rate[2] = mu;           // death from susceptible class
        rate[3] = gamma;        // recovery
        rate[4] = mu;           // death from infectious class
        rate[5] = mu;           // death from recovered class

        // compute the transition numbers
        trans[0] = rpois(rate[0]*dt);   // births are Poisson
        reulermultinom(2,S,&rate[1],dt,&trans[1]);
        reulermultinom(2,I,&rate[3],dt,&trans[3]);
        reulermultinom(1,R,&rate[5],dt,&trans[5]);

        // balance the equations
        S += trans[0]-trans[1]-trans[2];
        I += trans[1]-trans[3]-trans[4];
        R += trans[3]-trans[5];
      "),
      delta.t=1/52/20
      ),
      rinit=Csnippet("
        double m = pop/(S_0+I_0+R_0);
        S = nearbyint(m*S_0);
        I = nearbyint(m*I_0);
        R = nearbyint(m*R_0);
    "),
    paramnames=c("mu","pop","iota","gamma","Beta","sigma",
      "S_0","I_0","R_0"),
    statenames=c("S","I","R"),
    params=c(mu=1/50,iota=10,pop=50e6,gamma=26,Beta=400,sigma=0.1,
      S_0=0.07,I_0=0.001,R_0=0.93)
    ) -> ew1

  ew1 |>
    simulate() |>
    plot(variables=c("S","I","R"))

  ## A simple SIR model that tracks cumulative incidence.

  ew1 |>
    pomp(
      rprocess=euler(
        Csnippet("
        int nrate = 6;
        double rate[nrate];     // transition rates
        double trans[nrate];    // transition numbers
        double dW;

        // gamma noise, mean=dt, variance=(sigma^2 dt)
        dW = rgammawn(sigma,dt);

        // compute the transition rates
        rate[0] = mu*pop;       // birth into susceptible class
        rate[1] = (iota+Beta*I*dW/dt)/pop; // force of infection
        rate[2] = mu;           // death from susceptible class
        rate[3] = gamma;        // recovery
        rate[4] = mu;           // death from infectious class
        rate[5] = mu;           // death from recovered class

        // compute the transition numbers
        trans[0] = rpois(rate[0]*dt);   // births are Poisson
        reulermultinom(2,S,&rate[1],dt,&trans[1]);
        reulermultinom(2,I,&rate[3],dt,&trans[3]);
        reulermultinom(1,R,&rate[5],dt,&trans[5]);

        // balance the equations
        S += trans[0]-trans[1]-trans[2];
        I += trans[1]-trans[3]-trans[4];
        R += trans[3]-trans[5];
        H += trans[3];          // cumulative incidence
      "),
      delta.t=1/52/20
      ),
      rmeasure=Csnippet("
        double mean = H*rho;
        double size = 1/tau;
        reports = rnbinom_mu(size,mean);
    "),
    rinit=Csnippet("
        double m = pop/(S_0+I_0+R_0);
        S = nearbyint(m*S_0);
        I = nearbyint(m*I_0);
        R = nearbyint(m*R_0);
        H = 0;
    "),
    paramnames=c("mu","pop","iota","gamma","Beta","sigma","tau","rho",
      "S_0","I_0","R_0"),
    statenames=c("S","I","R","H"),
    params=c(mu=1/50,iota=10,pop=50e6,gamma=26,
      Beta=400,sigma=0.1,tau=0.001,rho=0.6,
      S_0=0.07,I_0=0.001,R_0=0.93)
    ) -> ew2

  ew2 |>
    simulate() |>
    plot()

  ## A simple SIR model that tracks weekly incidence.

  ew2 |>
    pomp(accumvars="H") -> ew3

  ew3 |>
    simulate() |>
    plot()

}
