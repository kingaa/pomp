setClass(
  "probed_pomp",
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
  definition=function (object, probes, params, nsim, seed = NULL, ...,
    verbose = getOption("verbose", FALSE))
  {

    if (missing(probes)) probes <- NULL
    if (missing(params)) params <- coef(object)
    if (missing(nsim)) nsim <- NULL

    probe.internal(
      object=object,
      probes=probes,
      params=params,
      nsim=nsim,
      seed=seed,
      ...,
      verbose=verbose
    )
  }
)

setMethod(
  "probe",
  signature=signature(object="probed_pomp"),
  definition=function (object, probes, nsim, seed = NULL, ...,
    verbose = getOption("verbose", FALSE)) {

    if (missing(probes)) probes <- object@probes
    if (missing(nsim)) nsim <- nrow(object@simvals)

    probe(
      as(object,"pomp"),
      probes=probes,
      nsim=nsim,
      seed=seed,
      ...,
      verbose=verbose
    )
  }
)

probe.internal <- function (object, probes, params, nsim, seed,
  .getnativesymbolinfo = TRUE, ..., verbose) {

  ep <- paste0("in ",sQuote("probe"),": ")
  verbose <- as.logical(verbose)

  object <- pomp(object,...)

  if (is.null(probes))
    stop(ep,sQuote("probes")," must be furnished.",call.=FALSE)
  if (!is.list(probes)) probes <- list(probes)
  if (!all(sapply(probes,is.function)))
    stop(ep,sQuote("probes")," must be a function or a list of functions.",
      call.=FALSE)
  if (!all(sapply(probes,function(f)length(formals(f))==1)))
    stop(ep,"each probe must be a function of a single argument.",call.=FALSE)

  if (is.list(params)) params <- unlist(params)
  if (is.null(params)) params <- numeric(0)
  if (length(params)==0)
    stop(ep,sQuote("params")," must be specified.",call.=FALSE)
  if (!is.numeric(params) || is.null(names(params)))
    stop(ep,sQuote("params")," must be furnished as a named numeric vector.",
      call.=FALSE)

  nsim <- as.integer(nsim)
  if (length(nsim) < 1)
    stop(ep,sQuote("nsim")," must be specified.",call.=FALSE)
  if (length(nsim) > 1 || !is.finite(nsim) || nsim <= 0)
    stop(ep,"number of simulations, ",sQuote("nsim"),
      ", must be a single positive integer.",call.=FALSE)

  seed <- as.integer(seed)

  gnsi <- as.logical(.getnativesymbolinfo)

  ## apply probes to data
  datval <- tryCatch(
    .Call(apply_probe_data,object,probes),
    error = function (e) {
      stop(ep,"applying probes to actual data: ",
        conditionMessage(e),call.=FALSE)
    }
  )

  nprobes <- length(datval)

  if (nprobes >= nsim)
    stop(ep,sQuote("nsim")," (=",nsim,"), should be (much) larger than the ",
      "number of probes (=",nprobes,").",call.=FALSE)

  pompLoad(object,verbose=verbose)
  on.exit(pompUnload(object,verbose=verbose))

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
      stop(ep,"applying probes to simulated data: ",
        conditionMessage(e),call.=FALSE)
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
      stop(ep,"in synthetic likelihood computation: ",  # nocov
        conditionMessage(e),call.=FALSE)                # nocov
    }
  )

  coef(object) <- params

  names(dimnames(simval)) <- c("rep","probe")

  new(
    "probed_pomp",
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
