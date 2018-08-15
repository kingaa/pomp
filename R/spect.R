## power spectrum
## Authors:
## Cai GoGwilt, Daniel Reuman, Aaron A. King

setClass(
  "spectd_pomp",
  contains="pomp",
  slots=c(
    kernel.width="numeric",
    transform.data="function",
    vars="character",
    freq="numeric",
    datspec="array",
    simspec="array",
    pvals="numeric",
    detrend="character"
  )
)

setGeneric("spect",function(object,...)standardGeneric("spect"))

setMethod(
  "spect",
  signature(object="pomp"),
  function (object, params, vars, kernel.width, nsim, seed = NULL,
    transform.data = identity, detrend = c("none","mean","linear","quadratic"),
    ...) {

    detrend <- match.arg(detrend)
    transform.data <- match.fun(transform.data)

    spect.internal(
      object,
      params=params,
      vars=vars,
      kernel.width=kernel.width,
      nsim=nsim,
      seed=seed,
      transform.data=transform.data,
      detrend=detrend,
      ...
    )

  }
)

setMethod(
  "spect",
  signature=signature(object="spectd_pomp"),
  definition=function (object, params, vars, kernel.width,
    nsim, seed = NULL, transform.data, detrend, ...) {

    if (missing(params)) params <- coef(object)
    if (missing(vars)) vars <- colnames(object@datspec)
    if (missing(kernel.width)) kernel.width <- object@kernel.width
    if (missing(nsim)) nsim <- nrow(object@simspec)
    if (missing(transform.data)) transform.data <- object@transform.data
    if (missing(detrend)) detrend <- object@detrend

    spect(
      as(object,"pomp"),
      params=params,
      vars=vars,
      kernel.width=kernel.width,
      nsim=nsim,
      seed=seed,
      transform.data=transform.data,
      detrend=detrend,
      ...
    )

  }
)

setMethod(
  "spect",
  signature=signature(object="missing"),
  definition=function (...) {
    stop("in ",sQuote("spect"),": ",sQuote("object")," is a required argument",call.=FALSE)
  }
)

setMethod(
  "spect",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("spect")," is not defined when ",sQuote("object")," is of class ",sQuote(class(object)),call.=FALSE)
  }
)

spect.internal <- function (object, params, vars, kernel.width, nsim,
  seed = NULL, transform.data, detrend, ...) {

  ep <- paste0("in ",sQuote("spect"),": ")

  if (missing(params)) params <- coef(object)
  if (is.list(params)) params <- unlist(params)

  if (missing(vars)) vars <- rownames(object@data)

  if (missing(kernel.width) || length(kernel.width) > 1 ||
      !is.numeric(kernel.width) ||
      !is.finite(kernel.width) || kernel.width < 0)
    stop(ep,sQuote("kernel.width"),
      " must be specified as a single positive integer.",call.=FALSE)

  if (missing(nsim) || length(nsim) > 1 || !is.numeric(nsim)||
      !is.finite(nsim) || (nsim<1))
    stop(ep,sQuote("nsim")," must be specified as a positive integer.",
      call.=FALSE)

  ker <- reuman.kernel(kernel.width)

  object <- pomp(object,...)

  pompLoad(object)
  on.exit(pompUnload(object))

  ds <- compute.spect.data(
    object,
    vars=vars,
    transform.data=transform.data,
    detrend=detrend,
    ker=ker
  )
  freq <- ds$freq
  datspec <- ds$spec

  simspec <- compute.spect.sim(
    object,
    params=params,
    vars=vars,
    nsim=nsim,
    seed=seed,
    transform.data=transform.data,
    detrend=detrend,
    ker=ker
  )

  pvals <- numeric(length(vars)+1)
  names(pvals) <- c(vars,"all")
  mean.simspec <- colMeans(simspec) # mean spectrum of simulations
  totdatdist <- 0
  totsimdist <- 0
  for (j in seq_along(vars)) {
    ## L-2 distance between data and mean simulated spectrum
    datdist <- sum((datspec[,j]-mean.simspec[,j])^2)
    ## L-2 distance betw. each sim. and mean simulated spectrum
    simdist <- vapply(
      seq_len(nsim),
      function(k)sum((simspec[k,,j]-mean.simspec[,j])^2),
      numeric(1)
    )
    pvals[j] <- (nsim+1-sum(simdist<datdist))/(nsim+1)
    totdatdist <- totdatdist+datdist
    totsimdist <- totsimdist+simdist
  }
  pvals[length(vars)+1] <- (nsim+1-sum(totsimdist<totdatdist))/(nsim+1)

  coef(object) <- params

  new(
    "spectd_pomp",
    object,
    kernel.width=kernel.width,
    transform.data=transform.data,
    vars=vars,
    detrend=detrend,
    freq=freq,
    datspec=datspec,
    simspec=simspec,
    pvals=pvals
  )
}

## detrends in one of several ways, according to type.
## tseries is a numeric vector,
pomp.detrend <- function (tseries, type) {
  switch(
    type,
    mean=tseries-mean(tseries),
    linear={
      m <- cbind(1,seq_along(tseries))
      .lm.fit(m,tseries)$residuals
    },
    quadratic={
      x <- seq_along(tseries)
      m <- cbind(1,x,x*x)
      .lm.fit(m,tseries)$residuals
    },
    tseries
  )
}

## The default smoothing kernel for the R spec.pgram function is weird.
## This function creates a better one.
reuman.kernel <- function (kernel.width) {
  ker <- kernel("modified.daniell",m=kernel.width)
  x <- seq.int(from=0,to=kernel.width,by=1)/kernel.width
  ker[[1L]] <- (15/(16*2*pi))*((x-1)^2)*((x+1)^2)
  ker[[1L]] <- ker[[1L]]/(2*sum(ker[[1L]][-1])+ker[[1L]][1L])
  attr(ker,"name") <- NULL
  ker
}

compute.spect.data <- function (object, vars, transform.data, detrend, ker) {

  ep <- paste0("in ",sQuote("spect"),": ")

  dat <- obs(object,vars)
  if (any(!is.finite(dat)))
    stop(ep,"missing or infinite values in the data.",call.=FALSE)

  dt <- diff(time(object,t0=FALSE))
  base.freq <- 1/mean(dt)
  dt.tol <- 0.025
  if (max(dt)-min(dt)>dt.tol*mean(dt))
    stop(ep,sQuote("spect")," assumes evenly spaced times.",call.=FALSE)

  for (j in seq_along(vars)) {
    sp <- spec.pgram(
      pomp.detrend(transform.data(dat[j,]),type=detrend),spans=ker,taper=0,
      pad=0,fast=FALSE,detrend=FALSE,plot=FALSE
    )
    if (j==1) {
      freq <- base.freq*sp$freq
      datspec <- array(
        dim=c(length(freq),nrow(dat)),
        dimnames=list(NULL,vars)
      )
    }
    datspec[,j] <- log10(sp$spec)
  }
  list(freq=freq,spec=datspec)
}

compute.spect.sim <- function (object, params, vars, nsim, seed, transform.data,
  detrend, ker) {

  ep <- paste0("in ",sQuote("spect"),": ")

  sims <- tryCatch(
    simulate(object,nsim=nsim,seed=seed,params=params,states=FALSE,obs=TRUE),
    error = function (e) {
      stop(ep,conditionMessage(e),call.=FALSE)
    }
  )

  sims <- sims[vars,,,drop=FALSE]

  if (any(!is.finite(sims)))
    stop(ep,"missing or infinite values in simulated data.",call.=FALSE)

  nobs <- length(vars)
  for (j in seq_len(nobs)) {
    for (k in seq_len(nsim)) {
      sp <- spec.pgram(pomp.detrend(transform.data(sims[j,k,]),type=detrend),
        spans=ker,taper=0,pad=0,fast=FALSE,detrend=FALSE,plot=FALSE)
      if ((j==1)&&(k==1)) {
        simspec <- array(
          dim=c(nsim,length(sp$freq),nobs),
          dimnames=list(NULL,NULL,vars)
        )
      }
      simspec[k,,j] <- log10(sp$spec)
    }
  }
  simspec
}
