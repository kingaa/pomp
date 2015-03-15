setClass(
         "spect.matched.pomp",
         contains="spect.pomp",
         slots=c(
           est="character",
           fail.value="numeric",
           weights="numeric",
           value="numeric",
           evals="integer",
           convergence="integer",
           msg="character"
           )
         )

spect.mismatch <- function (par, est, object, params,
                            vars, ker, nsim, seed,
                            transform, detrend, weights,
                            data.spec, fail.value) {
  if (missing(par)) par <- numeric(0)
  if (missing(est)) est <- integer(0)
  if (missing(params)) params <- coef(object)
  
  pompLoad(object)

  params[est] <- par
  
  ## vector of frequencies and estimated power spectum of data
  freq <- data.spec$freq
  datval <- data.spec$spec

  ## estimate power spectra of simulations
  simvals <- compute.spect.sim(
                               object,
                               vars=vars,
                               params=params,
                               nsim=nsim,
                               seed=seed,
                               transform=transform,
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

  pompUnload(object)
  mismatch
}

spect.match <- function(object, start, est = character(0),
                        vars, nsim, seed = NULL,
                        kernel.width, transform = identity, 
                        detrend = c("none","mean","linear","quadratic"),
                        weights,
                        method = c("subplex","Nelder-Mead","SANN"),
                        verbose = getOption("verbose"),
                        eval.only = FALSE, fail.value = NA, ...) {

  pompLoad(object)

  obj.fn <- spect.mismatch

  if (!is(object,"pomp"))
    stop(sQuote("object")," must be of class ",sQuote("pomp"))

  if (missing(start)) start <- coef(object)

  if (!eval.only&&(length(est)<1))
    stop("parameters to be estimated must be specified in ",sQuote("est"))
  if (!is.character(est)|!all(est%in%names(start)))
    stop(sQuote("est")," must refer to parameters named in ",sQuote("start"))
  par.index <- which(names(start)%in%est)
  
  if (missing(vars)) vars <- rownames(object@data)
            
  if (missing(nsim)) {
    if (is(object,"spect.pomp"))
      nsim <- nrow(object@simspec)
    else
      stop(sQuote("nsim")," must be supplied")
  }

  if (missing(kernel.width)) {
    if (is(object,"spect.pomp")) {
      kernel.width <- object@kernel.width
    } else {
      stop(sQuote("kernel.width")," must be specified")
    }
  }
  ker <- reuman.kernel(kernel.width)

  if (missing(transform)) {
    if (is(object,"spect.pomp")) {
      transform <- object@transform
    } else {
      transform <- identity
    }
  }

  if (nsim<1)
    stop(sQuote("nsim")," must be specified as a positive integer")

  if (missing(detrend)) {
    if (is(object,"spect.pomp")) {
      detrend <- object@detrend
    } else {
      detrend <- "none"
    }
  }
  detrend <- match.arg(detrend)

  method <- match.arg(method)
    
  ds <- compute.spect.data(
                           object,
                           vars=vars,
                           transform=transform,
                           detrend=detrend,
                           ker=ker
                           )

  if (missing(weights)) weights <- 1
  if (is.numeric(weights)) {
    if (length(weights)==1)
      weights <- rep(weights,length(ds$freq))
    if ((length(weights)!=length(ds$freq)))
      stop("if ",sQuote("weights")," is provided as a vector, it must have length ",length(ds$freq))
  } else if (is.function(weights)) {
    weights <- vapply(ds$freq,weights,numeric(1))
  } else {
    stop(sQuote("weights")," must be specified as a vector or as a function")
  }
  if (any((!is.finite(weights))|(weights<0)))
    stop(sQuote("weights")," should be nonnegative and finite")
  weights <- weights/mean(weights)

  fail.value <- as.numeric(fail.value)

  params <- start
  guess <- params[par.index]

  if (eval.only) {
    val <- obj.fn(
                  par=guess,
                  est=par.index,
                  object=object,
                  params=params,
                  vars=vars,
                  ker=ker,
                  nsim=nsim,
                  seed=seed,
                  transform=transform,
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
                              est=par.index,
                              object=object,
                              params=params,
                              vars=vars,
                              ker=ker,
                              nsim=nsim,
                              seed=seed,
                              transform=transform,
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
                   est=par.index,
                   object=object,
                   params=params,
                   vars=vars,
                   ker=ker,
                   nsim=nsim,
                   seed=seed,
                   transform=transform,
                   detrend=detrend,
                   weights=weights,
                   data.spec=ds,
                   fail.value=fail.value,
                   method=method,
                   control=list(...)
                   )
    }
    val <- opt$value
    params[par.index] <- opt$par
    conv <- opt$convergence
    evals <- opt$counts
    msg <- opt$message
  }

  pompUnload(object)

  new(
      "spect.matched.pomp",
      spect(
            object,
            params=params,
            vars=vars,
            kernel.width=kernel.width,
            nsim=nsim,
            seed=seed,
            transform=transform,
            detrend=detrend
            ),
      est=names(start)[par.index],
      fail.value=as.numeric(fail.value),
      value=val,
      weights=weights,
      convergence=as.integer(conv),
      evals=as.integer(evals),
      msg=as.character(msg)
      )
}

setMethod(
          "summary",
          "spect.matched.pomp",
          function (object, ...) {
            c(
              summary(as(object,"spect.pomp")),
              list(
                   est=object@est,
                   value=object@value,
                   eval=object@evals,
                   convergence=object@convergence
                   ),
              if(length(object@msg)>0) list(msg=object@msg) else NULL
              )
          }
          )
