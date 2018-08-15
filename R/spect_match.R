setClass(
  "spect_matched_pomp",
  contains="spectd_pomp",
  slots=c(
    est="character",
    fail.value="numeric",
    weights="numeric",
    method="character",
    value="numeric",
    evals="integer",
    convergence="integer",
    msg="character"
  )
)

setGeneric("spect.match",function(object,...)
  standardGeneric("spect.match"))

setMethod(
  "spect.match",
  signature=signature(object="pomp"),
  definition=function(object, start, est = character(0),
    vars, nsim, seed = NULL,
    kernel.width, transform.data,
    detrend = c("none","mean","linear","quadratic"),
    weights = 1,
    method = c("subplex","Nelder-Mead","SANN"),
    verbose = getOption("verbose"),
    fail.value = NA, ...) {
    ep <- paste0("in ",sQuote("spect.match"),": ")
    if (missing(start)) start <- coef(object)
    if (missing(vars)) vars <- rownames(object@data)
    if (missing(nsim)) stop(ep,sQuote("nsim")," must be supplied",call.=FALSE)
    if (missing(kernel.width)) stop(ep,sQuote("kernel.width")," must be specified",
      call.=FALSE)
    if (missing(transform.data)) transform.data <- identity
    transform.data <- match.fun(transform.data)
    detrend <- match.arg(detrend)
    method <- match.arg(method)

    spect.match.internal(object, start, est, vars, nsim, seed,
      kernel.width, transform.data, detrend, weights,
      method, verbose, fail.value, ...)
  }
)

setMethod(
  "spect.match",
  signature=signature(object="spectd_pomp"),
  definition=function(object, start, est = character(0), vars, nsim,
    seed = NULL, kernel.width, transform.data,
    detrend, weights = 1,
    method = c("subplex","Nelder-Mead","SANN"),
    verbose = getOption("verbose"), fail.value = NA,
    ...) {

    if (missing(start)) start <- object@params
    if (missing(vars)) vars <- object@vars
    if (missing(nsim)) nsim <- nrow(object@simspec)
    if (missing(kernel.width)) kernel.width <- object@kernel.width
    if (missing(transform.data)) transform.data <- object@transform.data
    if (missing(detrend)) detrend <- object@detrend
    method <- match.arg(method)

    spect.match.internal(object,start=start,est=est,vars=vars,nsim=nsim,
      seed=seed,kernel.width=kernel.width,
      transform.data=transform.data,detrend=detrend,
      weights=weights,method=method,
      verbose=verbose,fail.value=fail.value,...)

  }
)

setMethod(
  "spect.match",
  signature=signature(object="spect_matched_pomp"),
  definition=function(object, start, est, vars, nsim, seed = NULL,
    kernel.width, transform.data, detrend, weights, method,
    verbose = getOption("verbose"), fail.value,
    ...) {

    if (missing(start)) start <- object@params
    if (missing(est)) est <- object@est
    if (missing(vars)) vars <- object@vars
    if (missing(nsim)) nsim <- nrow(object@simspec)
    if (missing(kernel.width)) kernel.width <- object@kernel.width
    if (missing(transform.data)) transform.data <- object@transform.data
    if (missing(detrend)) detrend <- object@detrend
    if (missing(weights)) weights <- object@weights
    if (missing(method)) method <- object@method
    if (missing(fail.value)) fail.value <- object@fail.value

    spect.match.internal(object,start=start,est=est,vars=vars,nsim=nsim,
      seed=seed,kernel.width=kernel.width,
      transform.data=transform.data,detrend=detrend,
      weights=weights,method=method,
      verbose=verbose,fail.value=fail.value,...)

  }
)

spect.match.internal <- function(object, start, est, vars, nsim, seed = NULL,
  kernel.width, transform.data, detrend, weights,
  method, verbose, fail.value, ...) {

  ep <- paste0("in ",sQuote("spect.match"),": ")

  obj.fn <- spect.mismatch

  est <- as.character(est)
  eval.only <- (length(est)<1) || est=="" || is.na(est)
  if (!eval.only && !all(est %in% names(start)))
    stop(ep,sQuote("est")," must refer to parameters named in ",sQuote("start"),
      call.=FALSE)
  est.index <- which(names(start)%in%est)

  vars <- as.character(vars)
  if (!all(vars %in% rownames(object@data)))
    stop(ep,sQuote("vars")," must name data variables",call.=FALSE)

  nsim <- as.integer(nsim)
  if (length(nsim)<1 || !is.finite(nsim) || (nsim<1L))
    stop(ep,sQuote("nsim")," must be specified as a positive integer",call.=FALSE)

  ker <- reuman.kernel(kernel.width)

  pompLoad(object,verbose=verbose)
  on.exit(pompUnload(object,verbose=verbose))

  ds <- compute.spect.data(object,vars=vars,transform.data=transform.data,
    detrend=detrend,ker=ker)

  if (is.numeric(weights)) {
    if (length(weights)==1) weights <- rep(weights,length(ds$freq))
    if ((length(weights)!=length(ds$freq)))
      stop(ep,"if ",sQuote("weights")," is provided as a vector, it must have length ",
        length(ds$freq),call.=FALSE)
  } else if (is.function(weights)) {
    weights <- tryCatch(
      vapply(ds$freq,weights,numeric(1)),
      error = function (e) {
        stop(ep,"problem with ",sQuote("weights")," function: ",
          conditionMessage(e),call.=FALSE)
      }
    )
  } else {
    stop(ep,sQuote("weights"),
      " must be specified as a vector or as a function",call.=FALSE)
  }
  if (any((!is.finite(weights)) | (weights<0)))
    stop(ep,sQuote("weights")," should be nonnegative and finite",call.=FALSE)
  weights <- weights/mean(weights)

  fail.value <- as.numeric(fail.value)

  params <- start
  if (is.list(params)) params <- unlist(params)

  guess <- params[est.index]

  if (eval.only) {
    val <- obj.fn(
      par=guess,
      est=est.index,
      object=object,
      params=params,
      vars=vars,
      ker=ker,
      nsim=nsim,
      seed=seed,
      transform.data=transform.data,
      detrend=detrend,
      weights=weights,
      data.spec=ds,
      fail.value=fail.value
    )
    conv <- NA
    evals <- as.integer(c(1,0))
    msg <- "no optimization performed"
  } else {
    if (method == 'subplex') {
      opt <- subplex::subplex(
        par=guess,
        fn=obj.fn,
        est=est.index,
        object=object,
        params=params,
        vars=vars,
        ker=ker,
        nsim=nsim,
        seed=seed,
        transform.data=transform.data,
        detrend=detrend,
        weights=weights,
        data.spec=ds,
        fail.value=fail.value,
        control=list(...)
      )
    } else {
      opt <- optim(
        par=guess,
        fn=obj.fn,
        est=est.index,
        object=object,
        params=params,
        vars=vars,
        ker=ker,
        nsim=nsim,
        seed=seed,
        transform.data=transform.data,
        detrend=detrend,
        weights=weights,
        data.spec=ds,
        fail.value=fail.value,
        method=method,
        control=list(...)
      )
    }
    val <- opt$value
    params[est.index] <- opt$par
    conv <- opt$convergence
    evals <- opt$counts
    msg <- opt$message
  }

  new(
    "spect_matched_pomp",
    spect(
      object,
      params=params,
      vars=vars,
      kernel.width=kernel.width,
      nsim=nsim,
      seed=seed,
      transform.data=transform.data,
      detrend=detrend
    ),
    est=names(start)[est.index],
    vars=vars,
    fail.value=as.numeric(fail.value),
    value=val,
    weights=weights,
    method=method,
    convergence=as.integer(conv),
    evals=as.integer(evals),
    msg=as.character(msg)
  )
}

spect.mismatch <- function (par, est, object, params,
  vars, ker, nsim, seed,
  transform.data, detrend, weights,
  data.spec, fail.value) {

  params[est] <- par

  ## vector of frequencies and estimated power spectum of data
  freq <- data.spec$freq
  datval <- data.spec$spec

  pompLoad(object)
  on.exit(pompUnload(object))

  ## estimate power spectra of simulations
  simvals <- compute.spect.sim(
    object,
    vars=vars,
    params=params,
    nsim=nsim,
    seed=seed,
    transform.data=transform.data,
    detrend=detrend,
    ker=ker
  )
  ## simvals is an nsim x nfreq x nobs array

  ## compute a measure of the discrepancies between simulations and data
  discrep <- array(dim=c(length(freq),length(vars)))
  sim.means <- colMeans(simvals)
  for (j in seq_along(freq)) {
    for (k in seq_along(vars)) {
      discrep[j,k] <- ((datval[j,k]-sim.means[j,k])^2)/mean((simvals[,j,k]-sim.means[j,k])^2)
    }
    discrep[j,] <- weights[j]*discrep[j,]
  }

  if (!all(is.finite(discrep))) {
    mismatch <- fail.value
  } else {
    mismatch <- sum(discrep)
  }

  mismatch
}
