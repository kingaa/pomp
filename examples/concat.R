gompertz(sigma=2,tau=1) -> g
Np <- c(low=100,med=1000,high=10000)
lapply(
  Np,
  \(np) pfilter(g,Np=np)
) |>
  concat() -> pfs

pfs
coef(pfs)
logLik(pfs)
eff_sample_size(pfs)
cond_logLik(pfs)

pfs |> plot()
