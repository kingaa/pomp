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

setMethod(
  "summary",
  "spectd_pomp",
  function (object, ...) {
    list(
      coef=coef(object),
      nsim=nrow(object@simspec),
      pvals=object@pvals
    )
  }
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

setMethod(
  "plot",
  "spectd_pomp",
  function (x, max.plots.per.page = 4,
    plot.data = TRUE,
    quantiles = c(.025, .25, .5, .75, .975),
    quantile.styles = list(lwd=1, lty=1, col="gray70"),
    data.styles = list(lwd=2, lty=2, col="black")) {

    plot_spect.internal(x,max.plots.per.page=max.plots.per.page,
      plot.data=plot.data,quantiles=quantiles,quantile.styles=quantile.styles,
      data.styles=data.styles)

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

plot_spect.internal <- function (x, max.plots.per.page, plot.data,
  quantiles, quantile.styles, data.styles) {

  ep <- paste0("in ",sQuote("plot"),": ")

  spomp <- x
  nquants <- length(quantiles)

  if (!is.list(quantile.styles))
    stop(ep,sQuote("quantile.styles")," must be a list.",call.=FALSE)

  for (i in c("lwd", "lty", "col")) {
    if (is.null(quantile.styles[[i]]))
      quantile.styles[[i]] <- rep(1,nquants)
    if (length(quantile.styles[[i]])==1)
      quantile.styles[[i]] <- rep(quantile.styles[[i]],nquants)
    if (length(quantile.styles[[i]])<nquants) {
      warning(ep,sQuote("quantile.styles"),
        " contains an element with more than 1 entry but fewer entries than quantiles",
        call.=FALSE)
      quantile.styles[[i]]<-rep(quantile.styles[[i]],nquants)
    }
  }

  if (plot.data) {
    nreps <- ncol(spomp@datspec)

    if (!is.list(data.styles))
      stop(ep,sQuote("data.styles")," must be a list",call.=FALSE)

    for (i in c("lwd", "lty", "col")) {
      if(is.null(data.styles[[i]]))
        data.styles[[i]] <- rep(2,nreps)
      if(length(data.styles[[i]])==1)
        data.styles[[i]] <- rep(data.styles[[i]],nreps)
      if(length(data.styles[[i]]) < nreps) {
        warning(ep,sQuote("data.styles"),
          "contains an element with more than 1 entry but ",
          "fewer entries than observed variables",
          call.=FALSE)
        data.styles[[i]] <- rep(data.styles[[i]],nreps)
      }
    }
  }

  dimsim <- dim(spomp@simspec)
  nfreq <- dimsim[2]
  nobs <- dimsim[3]
  oldpar <- par(
    mfrow=c(min(nobs,max.plots.per.page),1),
    mar=c(3,3,1,0.5),
    mgp=c(2,1,0),
    bty="l",
    ask=if (nobs>max.plots.per.page) TRUE else par("ask")
  )
  on.exit(par(oldpar))
  ylabs <- dimnames(spomp@simspec)[[3]]
  for (i in seq_len(nobs)) {
    spectraquants <- array(dim=c(nfreq,length(quantiles)))
    for (j in seq_len(nfreq))
      spectraquants[j,] <- quantile(
        spomp@simspec[,j,i],
        probs=quantiles
      )
    if (plot.data) {
      ylimits <- c(
        min(spectraquants,spomp@datspec[,i]),
        max(spectraquants,spomp@datspec[,i])
      )
    } else {
      ylimits <- c(
        min(spectraquants),
        max(spectraquants)
      )
    }
    plot(
      NULL,
      xlim=range(spomp@freq),ylim=ylimits,
      xlab=if (i==nobs) "frequency" else "",
      ylab=expression(paste(log[10],"power"))
    )
    title(
      main=paste0(ylabs[i],", p = ",round(spomp@pvals[i],4)),
      line=0
    )
    for (j in seq_along(quantiles)) {
      lines(
        x=spomp@freq,
        y=spectraquants[,j],
        lwd=quantile.styles$lwd[j],
        lty=quantile.styles$lty[j],
        col=quantile.styles$col[j]
      )
    }
    if(plot.data) {
      lines(
        x=spomp@freq,
        y=spomp@datspec[,i],
        lty=data.styles$lty[i],
        lwd=data.styles$lwd[i],
        col=data.styles$col[i]
      )
    }
  }
}
