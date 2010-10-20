setClass(
         "probe.matched.pomp",
         contains="probed.pomp",
         representation=representation(
           weights="numeric",
           fail.value="numeric",
           evals="integer",
           value="numeric",
           convergence="integer",
           msg="character"
           )
         )

setMethod(
          "summary",
          "probe.matched.pomp",
          function (object, ...) {
            c(
              summary(as(object,"probed.pomp")),
              list(
                   weights=object@weights,
                   value=object@value,
                   eval=object@evals,
                   convergence=object@convergence
                   ),
              if(length(object@msg)>0) list(msg=object@msg) else NULL
              )
          }
          )

probe.mismatch <- function (par, est, object, probes, params,
                            nsim = 1, seed = NULL,
                            weights, datval,
                            fail.value = NA) {
  if (missing(par)) par <- numeric(0)
  if (missing(est)) est <- integer(0)
  if (missing(params)) params <- coef(object)
  
  params[est] <- par
  
  ## apply probes to model simulations
  simval <- .Call(
                  apply_probe_sim,
                  object=object,
                  nsim=nsim,
                  params=params,
                  seed=seed,
                  probes=probes,
                  datval=datval
                  )
  
  ## compute a measure of the discrepancies between simulations and data
  sim.means <- colMeans(simval)
  simval <- sweep(simval,2,sim.means)
  discrep <- ((datval-sim.means)^2)/colMeans(simval^2)
  if ((length(weights)>1) && (length(weights)!=length(discrep)))
    stop(length(discrep)," probes have been computed, but ",length(weights)," have been supplied")
  if (!all(is.finite(discrep))) {
    mismatch <- fail.value 
  } else if (length(weights)>1) {
    mismatch <- sum(discrep*weights)/sum(weights)
  } else {
    mismatch <- sum(discrep)
  }

  mismatch
}

neg.synth.loglik <- function (par, est, object, probes, params,
                              nsim = 1, seed = NULL,
                              weights, datval,
                              fail.value = NA) {
  if (missing(par)) par <- numeric(0)
  if (missing(est)) est <- integer(0)
  if (missing(params)) params <- coef(object)
  
  params[est] <- par
  
  ## apply probes to model simulations
  simval <- .Call(
                  apply_probe_sim,
                  object=object,
                  nsim=nsim,
                  params=params,
                  seed=seed,
                  probes=probes,
                  datval=datval
                  )
  
  ll <- .Call(synth_loglik,simval,datval)
  -ll
}

probe.match <- function(object, start, est = character(0),
                        probes, weights,
                        nsim, seed = NULL,
                        method = c("subplex","Nelder-Mead","SANN"),
                        verbose = getOption("verbose"), 
                        eval.only = FALSE, fail.value = NA, ...) {

  obj.fn <- probe.mismatch

  if (!is(object,"pomp"))
    stop(sQuote("object")," must be of class ",sQuote("pomp"))

  if (missing(start))
    start <- coef(object)

  if (!eval.only&&(length(est)<1))
    stop("parameters to be estimated must be specified in ",sQuote("est"))
  if (!is.character(est)|!all(est%in%names(start)))
    stop(sQuote("est")," must refer to parameters named in ",sQuote("start"))
  par.index <- which(names(start)%in%est)
  
  if (missing(probes)) {
    if (is(object,"probed.pomp"))
      probes <- object@probes
    else
      stop(sQuote("probes")," must be supplied")
  }
  if (!is.list(probes)) probes <- list(probes)
  if (!all(sapply(probes,is.function)))
    stop(sQuote("probes")," must be a function or a list of functions")
  if (!all(sapply(probes,function(f)length(formals(f))==1)))
    stop("each probe must be a function of a single argument")            

  if (missing(weights)) weights <- 1

  method <- match.arg(method)

  params <- start
  guess <- params[par.index]

  datval <- .Call(apply_probe_data,object,probes) # apply probes to data
  
  if (eval.only) {
    val <- obj.fun(
                   par=guess,
                   est=par.index,
                   object=object,
                   probes=probes,
                   params=params,
                   nsim=nsim,
                   seed=seed,
                   weights=weights,
                   datval=datval,
                   fail.value=fail.value
                   )
    conv <- 0
    evals <- as.integer(c(1,0))
    msg <- paste("no optimization performed")
  } else {
    if (method == 'subplex') {
      opt <- subplex::subplex(
                              par=guess,
                              fn=obj.fn,
                              est=par.index,
                              object=object,
                              probes=probes,
                              params=params,
                              nsim=nsim,
                              seed=seed,
                              weights=weights,
                              datval=datval,
                              fail.value=fail.value,
                              control=list(...)
                              )
    } else {
      opt <- optim(
                   par=guess,
                   fn=obj.fn,
                   est=par.index,
                   object=object,
                   probes=probes,
                   params=params,
                   nsim=nsim,
                   seed=seed,
                   weights=weights,
                   datval=datval,
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

  new(
      "probe.matched.pomp",
      probe(
            as(object,"pomp"),
            probes=probes,
            params=params,
            nsim=nsim,
            seed=seed
            ),
      weights=weights,
      fail.value=as.numeric(fail.value),
      value=val,
      convergence=as.integer(conv),
      evals=as.integer(evals),
      msg=as.character(msg)
      )
}
