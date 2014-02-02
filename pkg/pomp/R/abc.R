## define the abc class
setClass(
         'abc',
         contains='pomp',
         slots=c(
           pars = 'character',
           transform = 'logical',
           Nabc = 'integer',
           dprior = 'function',
           probes='list',
           scale = 'numeric',
           epsilon = 'numeric',
           hyperparams = 'list',
           random.walk.sd = 'numeric',
           conv.rec = 'matrix'
           )
         )

## ABC algorithm functions
setGeneric('abc',function(object,...)standardGeneric("abc"))

abc.internal <- function (object, Nabc,
                          start, pars, dprior.fun,
                          rw.sd, probes,
                          hyperparams,
                          epsilon, scale,
                          verbose, transform,
                          .ndone, .getnativesymbolinfo = TRUE,
                          ...) {

  gnsi <- as.logical(.getnativesymbolinfo)

  transform <- as.logical(transform)

  Nabc <- as.integer(Nabc)

  if (length(start)==0)
    stop(
         "abc error: ",sQuote("start")," must be specified if ",
         sQuote("coef(object)")," is NULL",
         call.=FALSE
         )

  start.names <- names(start)
  if (is.null(start.names))
    stop("abc error: ",sQuote("start")," must be a named vector",call.=FALSE)

  if (missing(rw.sd))
    stop("abc error: ",sQuote("rw.sd")," must be specified",call.=FALSE)
  rw.names <- names(rw.sd)
  if (is.null(rw.names) || any(rw.sd<0))
    stop("abc error: ",sQuote("rw.sd")," must be a named non-negative numerical vector",call.=FALSE)
  if (!all(rw.names%in%start.names))
    stop("abc error: all the names of ",sQuote("rw.sd")," must be names of ",sQuote("start"),call.=FALSE)
  rw.names <- names(rw.sd[rw.sd>0])
  if (length(rw.names) == 0)
    stop("abc error: ",sQuote("rw.sd")," must have one positive entry for each parameter to be estimated",call.=FALSE)

  if (
      !is.character(pars) ||
      !all(pars%in%start.names) ||
      !all(pars%in%rw.names)
      )
    stop(
         "abc error: ",
         sQuote("pars"),
         " must be a mutually disjoint subset of ",
         sQuote("names(start)"),
         " and must correspond to positive random-walk SDs specified in ",
         sQuote("rw.sd"),
         call.=FALSE
         )

  if (!all(rw.names%in%pars)) {
    extra.rws <- rw.names[!(rw.names%in%pars)]
    warning(
            "abc warning: the variable(s) ",
            paste(extra.rws,collapse=", "),
            " have positive random-walk SDs specified, but are not included in ",
            sQuote("pars"),
            ". These random walk SDs are ignored.",
            call.=FALSE
            )
  }
  rw.sd <- rw.sd[pars]
  rw.names <- names(rw.sd)

  if (!is.list(probes)) probes <- list(probes)
  if (!all(sapply(probes,is.function)))
    stop(sQuote("probes")," must be a function or a list of functions")
  if (!all(sapply(probes,function(f)length(formals(f))==1)))
    stop("each probe must be a function of a single argument")

  ntimes <- length(time(object))
  
  if (verbose) {
    cat("performing",Nabc,"ABC iteration(s) to estimate parameter(s)",
        paste(pars,collapse=", "))
    cat(" using random-walk with SD\n")
    print(rw.sd)
  }

  theta <- start
  log.prior <- dprior.fun(params=theta,hyperparams=hyperparams,log=TRUE)
  ## we suppose that theta is a "match", which does the right thing for continue() and
  ## should have negligible effect unless doing many short calls to continue()

  conv.rec <- matrix(
                     data=NA,
                     nrow=Nabc+1,
                     ncol=length(theta),
                     dimnames=list(
                       rownames=seq(from=0,to=Nabc,by=1),
                       colnames=names(theta)
                       )
                     )

  if (!all(is.finite(theta[pars]))) {
    stop(
         sQuote("abc"),
         " error: cannot estimate non-finite parameters: ",
         paste(
               pars[!is.finite(theta[pars])],
               collapse=","
               ),
         call.=FALSE
         )
  }

  po <- as(object,"pomp")
  
  ## apply probes to data
  datval <- try(.Call(apply_probe_data,po,probes),silent=FALSE)
  if (inherits(datval,'try-error'))
    stop("abc error: error in ",sQuote("apply_probe_data"),call.=FALSE)

  conv.rec[1,names(theta)] <- theta

  for (n in seq_len(Nabc)) { # main loop

    theta.prop <- theta
    theta.prop[pars] <- rnorm(n=length(pars),mean=theta.prop[pars],sd=rw.sd)

    ## compute the probes for the proposed new parameter values

    simval <- try(
                  .Call(
                        apply_probe_sim,
                        object=po,
                        nsim=1,
                        params=theta.prop,
                        seed=NULL,
                        probes=probes,
                        datval=datval
                        ),
                  silent=FALSE
                  )

    if (inherits(simval,'try-error'))
      stop("abc error: error in ",sQuote("apply_probe_sim"),call.=FALSE)

    ## ABC update rule
    distance <- sum(((datval-simval)/scale)^2)
    if( (is.finite(distance)) && (distance<epsilon) ){ 
      log.prior.prop <- dprior.fun(params=theta.prop,hyperparams=hyperparams,log=TRUE)
      if (runif(1) < exp(log.prior.prop-log.prior)) {
        theta <- theta.prop
        log.prior <- log.prior.prop
      }
    }

    ## store a record of this iteration
    conv.rec[n+1,names(theta)] <- theta
    if (verbose && (n%%5==0))
      cat("ABC iteration ",n," of ",Nabc," completed\n")

  }

  new(
      'abc',
      po,
      params=theta,
      pars = pars,
      transform = transform,
      Nabc = Nabc,
      dprior = dprior.fun,
      probes=probes,
      scale = scale,
      epsilon = epsilon,
      hyperparams = hyperparams,
      random.walk.sd = rw.sd,
      conv.rec=conv.rec
      )

}

setMethod(
          "abc",
          signature=signature(object="pomp"),
          function (object, Nabc = 1,
                    start, pars, rw.sd,
                    dprior, probes, scale, epsilon, hyperparams,
                    verbose = getOption("verbose"),
                    transform = FALSE,
                    ...) {
            
            transform <- as.logical(transform)

            if (missing(start)) start <- coef(object,transform=transform)

            if (missing(rw.sd))
              stop("abc error: ",sQuote("rw.sd")," must be specified",call.=FALSE)

            if (missing(pars)) {
              pars <- names(rw.sd)[rw.sd>0]
            }

            if (missing(probes))
              stop("abc error: ",sQuote("probes")," must be specified",call.=FALSE)

            if (missing(scale))
              stop("abc error: ",sQuote("scale")," must be specified",call.=FALSE)

            if (missing(epsilon))
              stop("abc error: ",sQuote("abc match criterion, epsilon,")," must be specified",call.=FALSE)

            if (missing(hyperparams))
              stop("abc error: ",sQuote("hyperparams")," must be specified",call.=FALSE)

            if (missing(dprior)) {         # use default flat improper prior
              dprior <- function (params, hyperparams, log) {
                if (log) 0 else 1
              }
            } else {
              dprior <- match.fun(dprior)
              if (!all(c('params','hyperparams','log')%in%names(formals(dprior))))
                stop(
                     "abc error: ",
                     sQuote("dprior"),
                     " must be a function of prototype ",
                     sQuote("dprior(params,hyperparams,log)"),
                     call.=FALSE
                     )
            }
            
            abc.internal(
                         object=object,
                         Nabc=Nabc,
                         start=start,
                         pars=pars,
                         dprior.fun=dprior,
                         probes=probes,
                         scale=scale,
                         epsilon=epsilon,
                         rw.sd=rw.sd,
                         hyperparams=hyperparams,
                         verbose=verbose,
                         transform=transform,
                         .ndone=0
                         )
          }
          )

setMethod(
          "abc",
          signature=signature(object="probed.pomp"),
          function (object, Nabc = 1,
                    start, pars, rw.sd,
                    dprior, probes, scale, epsilon, hyperparams,
                    verbose = getOption("verbose"),
                    transform = FALSE,
                    ...) {

            transform <- as.logical(transform)
            
            if (missing(start)) start <- coef(object,transform=transform)

            if (missing(rw.sd))
              stop("abc error: ",sQuote("rw.sd")," must be specified",call.=FALSE)

            if (missing(pars)) {
              pars <- names(rw.sd)[rw.sd>0]
            }

            if (missing(probes)) probes <- object@probes
            if (missing(scale)) probes <- object@scale
            if (missing(epsilon)) probes <- object@epsilon

            if (missing(hyperparams))
              stop("abc error: ",sQuote("hyperparams")," must be specified",call.=FALSE)

            if (missing(dprior)) {         # use default flat improper prior
              dprior <- function (params, hyperparams, log) {
                if (log) 0 else 1
              }
            } else {
              dprior <- match.fun(dprior)
              if (!all(c('params','hyperparams','log')%in%names(formals(dprior))))
                stop(
                     "abc error: ",
                     sQuote("dprior"),
                     " must be a function of prototype ",
                     sQuote("dprior(params,hyperparams,log)"),
                     call.=FALSE
                     )
            }
            
            abc.internal(
                         object=as(object,"pomp"),
                         Nabc=Nabc,
                         start=start,
                         pars=pars,
                         dprior.fun=dprior,
                         probes=probes,
                         scale=scale,
                         epsilon=epsilon,
                         rw.sd=rw.sd,
                         hyperparams=hyperparams,
                         verbose=verbose,
                         transform=transform,
                         .ndone=0
                         )
          }
          )

setMethod(
          "abc",
          signature=signature(object="abc"),
          function (object, Nabc,
                    start, pars, rw.sd,
                    dprior, probes, scale, epsilon, hyperparams,
                    verbose = getOption("verbose"),
                    transform,
                    ...) {

            if (missing(Nabc)) Nabc <- object@Nabc
            if (missing(start)) start <- coef(object)
            if (missing(pars)) pars <- object@pars
            if (missing(rw.sd)) pars <- object@random.walk.sd
            if (missing(dprior)) dprior <- object@dprior
            if (missing(probes)) probes <- object@probes
            if (missing(scale)) probes <- object@scale
            if (missing(epsilon)) probes <- object@epsilon
            if (missing(hyperparams)) hyperparams <- object@hyperparams
            if (missing(transform)) transform <- object@transform
            transform <- as.logical(transform)

            abc.internal(
                         object=as(object,"pomp"),
                         Nabc=Nabc,
                         start=start,
                         pars=pars,
                         dprior.fun=dprior,
                         rw.sd=rw.sd,
                         probes=probes,
                         scale=scale,
                         epsilon=epsilon,
                         hyperparams=hyperparams,
                         verbose=verbose,
                         transform=transform,
                         .ndone=0
                         )
          }
          )

setMethod(
          'continue',
          signature=signature(object='abc'),
          function (object, Nabc = 1, ...) {

            ndone <- object@Nabc
            
            obj <- abc(
                       object=object,
                       Nabc=Nabc,
                       .ndone=ndone,
                       ...
                       )
            
            obj@conv.rec <- rbind(
                                  object@conv.rec[,colnames(obj@conv.rec)],
                                  obj@conv.rec[-1,]
                                  )
            obj@Nabc <- as.integer(ndone+Nabc)
            obj
          }
          )
