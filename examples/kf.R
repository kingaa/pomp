\donttest{ # takes too long for R CMD check

  if (require(dplyr)) {

    gompertz() -> po

    po |>
      as.data.frame() |>
      mutate(
        logY=log(Y)
      ) |>
      select(time,logY) |>
      pomp(times="time",t0=0) |>
      kalmanFilter(
        X0=c(logX=0),
        A=matrix(exp(-0.1),1,1),
        Q=matrix(0.01,1,1),
        C=matrix(1,1,1),
        R=matrix(0.01,1,1)
      ) -> kf

    po |>
      pfilter(Np=1000) -> pf

    kf$logLik
    logLik(pf) + sum(log(obs(pf)))
    
  }
}
