library(pomp)
set.seed(964484044)
sir_step <- Csnippet("
  double dN_SI = rbinom(S,1-exp(-exp(Beta)*(I/10000)*dt));
  double dN_IR = rbinom(I,1-exp(-exp(mu_IR)*dt));
  S -= dN_SI;
  I += dN_SI - dN_IR;
  R += dN_IR;
  H += dN_SI;")
sir_init <- Csnippet("
  S = 10000-1;
  I = 1;
  R = 0;
  H = 0;")
dmeas <- Csnippet("
  lik = dbinom(reports,H,1 - exp(-exp(rho)),give_log);")
rmeas <- Csnippet("
  reports = rbinom(H,1-exp(-exp(rho)));")
simulate(
  times=1:10,
  t0=1,
  rmeasure = rmeas,
  dmeasure = dmeas,
  rprocess = discrete_time(step.fun = Csnippet(sir_step), delta.t = 1),
  obsnames="reports",
  statenames = c("S", "I", "R", "H"),
  paramnames = c("Beta", "mu_IR", "rho"),
  rinit = sir_init,
  params = c(Beta = log(1), mu_IR = log(0.1), rho = log(0.1))
) |>
  pmcmc(
    Np = 1000,
    Nmcmc = 100,
    proposal = mvn_diag_rw(c(Beta = 5, mu_IR = 5, rho = 5))
#    proposal = mvn_diag_rw(c(Beta = 0.2, mu_IR = 0.2, rho = 0.2))
  ) -> pmcmcSIR
