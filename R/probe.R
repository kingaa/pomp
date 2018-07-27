setClass(
  "probed.pomp",
  contains="pomp",
  slots=c(
    probes="list",
    datvals="numeric",
    simvals="array",
    quantiles="numeric",
    pvals="numeric",
    synth.loglik="numeric",
    seed="integer"
  )
)

setMethod(
  "probe",
  signature=signature(object="pomp"),
  definition=function (object, probes, params, nsim = 1, seed = NULL, ...)
  {
    probe.internal(
      object=object,
      probes=probes,
      params=params,
      nsim=nsim,
      seed=seed,
      ...
    )
  }
)

setMethod(
  "probe",
  signature=signature(object="probed.pomp"),
  definition=function (object, probes, nsim, seed, ...) {

    if (missing(probes)) probes <- object@probes
    if (missing(nsim)) nsim <- nrow(object@simvals)
    if (missing(seed)) seed <- object@seed

    probe(
      as(object,"pomp"),
      probes=probes,
      nsim=nsim,
      seed=seed,
      ...
    )
  }
)

probe.internal <- function (object, probes, params, nsim = 1L, seed = NULL,
  .getnativesymbolinfo = TRUE, ...) {

  ep <- paste0("in ",sQuote("probe"),": ")

  if (!is.list(probes)) probes <- list(probes)
  if (!all(sapply(probes,is.function)))
    stop(ep,sQuote("probes")," must be a function or a list of functions",call.=FALSE)
  if (!all(sapply(probes,function(f)length(formals(f))==1)))
    stop(ep,"each probe must be a function of a single argument",call.=FALSE)

  gnsi <- as.logical(.getnativesymbolinfo)

  seed <- as.integer(seed)

  if (missing(params)) params <- coef(object)
  if (is.list(params)) params <- unlist(params)

  pompLoad(object)

  ## apply probes to data
  datval <- tryCatch(
    .Call(apply_probe_data,object,probes),
    error = function (e) {
      stop(ep,"applying probes to actual data: ",conditionMessage(e),call.=FALSE)
    }
  )
  nprobes <- length(datval)

  ## apply probes to model simulations
  simval <- tryCatch(
    .Call(
      apply_probe_sim,
      object=object,
      nsim=nsim,
      params=params,
      seed=seed,
      probes=probes,
      datval=datval,
      gnsi=gnsi
    ),
    error = function (e) {
      stop(ep,"applying probes to simulated data: ",conditionMessage(e),call.=FALSE)
    }
  )

  pvals <- numeric(nprobes)
  names(pvals) <- names(datval)
  quants <- numeric(nprobes)
  names(quants) <- names(datval)
  for (k in seq_len(nprobes)) {
    r <- min(sum(simval[,k]>datval[k]),sum(simval[,k]<datval[k]))
    tails <- (r+1)/(nsim+1)
    pvals[k] <- min(2*tails,1)
    quants[k] <- sum(simval[,k]<datval[k])/nsim
  }

  ll <- tryCatch(
    .Call(synth_loglik,simval,datval),
    error = function (e) {
      stop(ep,"in synthetic likelihood computation: ",conditionMessage(e),call.=FALSE)
    }
  )

  coef(object) <- params

  pompUnload(object)

  new(
    "probed.pomp",
    object,
    probes=probes,
    datvals=datval,
    simvals=simval,
    quantiles=quants,
    pvals=pvals,
    synth.loglik=ll,
    seed=seed
  )
}
