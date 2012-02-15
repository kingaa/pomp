setClass(
         "probe.matched.pomp",
         contains="probed.pomp",
         representation=representation(
           est="character",
           weights="numeric",
           fail.value="numeric",
           value="numeric",
           evals="integer",
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
                   est=object@est,
                   value=object@value,
                   eval=object@evals,
                   convergence=object@convergence
                   ),
              if(length(object@msg)>0) list(msg=object@msg) else NULL
              )
          }
          )

probe.match.objfun <- function (object, params, est, probes,
                                nsim = 1, seed = NULL, fail.value = NA, ...) {

  if (missing(est)) est <- character(0)
  if (!is.character(est)) stop(sQuote("est")," must be a vector of parameter names")
  if (missing(params)) params <- coef(object)
  if ((!is.numeric(params))||(is.null(names(params))))
    stop(sQuote("params")," must be a named numeric vector")
  par.est.idx <- match(est,names(params))
  if (any(is.na(par.est.idx)))
    stop("parameter(s): ",sQuote(est[is.na(par.est.idx)])," not found in ",sQuote("params"))
  
  if (!is.list(probes)) probes <- list(probes)
  if (!all(sapply(probes,is.function)))
    stop(sQuote("probes")," must be a function or a list of functions")
  if (!all(sapply(probes,function(f)length(formals(f))==1)))
    stop("each probe must be a function of a single argument")            

  datval <- .Call(apply_probe_data,object,probes) # apply probes to data
  nprobes <- length(datval)
  if (nprobes > nsim)
    stop(sQuote("nsim"),"(=",nsim,") should be (much) larger than the number of probes (=",nprobes,")")
    
  obj.fun <- function (par) {
    params[par.est.idx] <- par
  
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
    if (is.finite(ll)||is.na(fail.value)) -ll else fail.value
  }

  obj.fun
}

probe.match.internal <- function(object, start, est,
                                 probes, weights,
                                 nsim, seed,
                                 method, verbose,
                                 eval.only, fail.value, ...) {

  if (eval.only) {
    est <- character(0)
    guess <- numeric(0)
  } else {
    if (!is.character(est)) stop(sQuote("est")," must be a vector of parameter names")
    if (length(start)<1)
      stop(sQuote("start")," must be supplied if ",sQuote("object")," contains no parameters")
    if (is.null(names(start))||(!all(est%in%names(start))))
      stop(sQuote("est")," must refer to parameters named in ",sQuote("start"))
    guess <- start[est]
  }

  obj <- as(object,"pomp")
  coef(obj) <- start

  obj.fn <- probe.match.objfun(
                               obj,
                               est=est,
                               probes=probes,
                               nsim=nsim,
                               seed=seed,
                               fail.value=fail.value
                               )

  
  if (eval.only) {

    val <- obj.fn(guess)
    conv <- NA
    evals <- as.integer(c(1,0))
    msg <- paste("no optimization performed")

  } else {

    if (method == 'subplex') {
      opt <- subplex::subplex(par=guess,fn=obj.fn,control=list(...))
    } else if (method=="sannbox") {
      opt <- sannbox(par=guess,fn=obj.fn,control=list(...))
    } else {
      opt <- optim(par=guess,fn=obj.fn,method=method,control=list(...))
    }

    if (!is.null(names(opt$par)) && !all(est==names(opt$par)))
      stop("mismatch between parameter names returned by optimizer and ",sQuote("est"))
    coef(obj,est) <- unname(opt$par)
    msg <- if (is.null(opt$message)) character(0) else opt$message
    conv <- opt$convergence
    val <- opt$value
    evals <- opt$counts
  }

  new(
      "probe.matched.pomp",
      probe(
            obj,
            probes=probes,
            nsim=nsim,
            seed=seed
            ),
      est=as.character(est),
      value=val,
      convergence=as.integer(conv),
      evals=as.integer(evals),
      msg=as.character(msg)
      )
}

setGeneric("probe.match",function(object,...)standardGeneric("probe.match"))

setMethod(
          "probe.match",
          signature=signature(object="pomp"),
          function(object, start, est = character(0),
                   probes, weights,
                   nsim, seed = NULL,
                   method = c("subplex","Nelder-Mead","SANN","BFGS","sannbox"),
                   verbose = getOption("verbose"), 
                   eval.only = FALSE, fail.value = NA, ...) {
            
            if (missing(start)) start <- coef(object)

            if (missing(probes))
              stop(sQuote("probes")," must be supplied")

            if (missing(nsim))
              stop(sQuote("nsim")," must be supplied")

            if (missing(weights)) weights <- 1

            method <- match.arg(method)
            
            probe.match.internal(
                                 object=object,
                                 start=start,
                                 est=est,
                                 probes=probes,
                                 weights=weights,
                                 nsim=nsim,
                                 seed=seed,
                                 method=method,
                                 verbose=verbose,
                                 eval.only=eval.only,
                                 fail.value=fail.value,
                                 ...
                                 )
          }
          )

setMethod(
          "probe.match",
          signature=signature(object="probed.pomp"),
          function(object, start, est = character(0),
                   probes, weights,
                   nsim, seed = NULL,
                   method = c("subplex","Nelder-Mead","SANN","BFGS","sannbox"),
                   verbose = getOption("verbose"), 
                   eval.only = FALSE, fail.value = NA, ...) {
            
            if (missing(start)) start <- coef(object)

            if (missing(probes))
              probes <- object@probes

            if (missing(nsim))
              nsim <- nrow(object@simvals)
            
            if (missing(weights)) weights <- 1

            method <- match.arg(method)
            
            probe.match.internal(
                                 object=object,
                                 start=start,
                                 est=est,
                                 probes=probes,
                                 weights=weights,
                                 nsim=nsim,
                                 seed=seed,
                                 method=method,
                                 verbose=verbose,
                                 eval.only=eval.only,
                                 fail.value=fail.value,
                                 ...
                                 )
          }
          )

setMethod(
          "probe.match",
          signature=signature(object="probe.matched.pomp"),
          function(object, start, est,
                   probes, weights,
                   nsim, seed = NULL,
                   method = c("subplex","Nelder-Mead","SANN","BFGS","sannbox"),
                   verbose = getOption("verbose"), 
                   eval.only = FALSE, fail.value, ...) {
            
            if (missing(start)) start <- coef(object)

            if (missing(est)) est <- object@est

            if (missing(probes))
              probes <- object@probes

            if (missing(nsim))
              nsim <- nrow(object@simvals)
            
            if (missing(weights)) weights <- 1

            if (missing(fail.value)) fail.value <- object@fail.value

            method <- match.arg(method)
            
            probe.match.internal(
                                 object=object,
                                 start=start,
                                 est=est,
                                 probes=probes,
                                 weights=weights,
                                 nsim=nsim,
                                 seed=seed,
                                 method=method,
                                 verbose=verbose,
                                 eval.only=eval.only,
                                 fail.value=fail.value,
                                 ...
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
  if (is.finite(ll)||is.na(fail.value)) -ll else fail.value
}
