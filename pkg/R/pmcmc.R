## define the pmcmc class
setClass(
         'pmcmc',
         representation(
                        pars = 'character',
                        Nmcmc = 'integer',
                        dprior = 'function',
                        hyperparams = 'list',
                        Np = 'integer',
                        random.walk.sd = 'numeric',
                        filter.mean = 'matrix',
                        conv.rec = 'matrix',
                        eff.sample.size = 'numeric',
                        cond.loglik = 'numeric',
                        loglik = 'numeric',
                        log.prior = 'numeric'
                        ),
         contains='pomp'
         )

## PMCMC algorithm functions
pmcmc <- function (object, ... )
  stop("function ",sQuote("pmcmc")," is undefined for objects of class ",sQuote(class(object)))
setGeneric('pmcmc')

dprior <- function (object, params, log)
  stop("function ",sQuote("dprior")," is undefined for objects of class ",sQuote(class(object)))
setGeneric('dprior')

setMethod(
          "dprior",
          "pmcmc",
          function (object, params, log = FALSE) {
            do.call(object@dprior,list(params=params,hyperparams=object@hyperparams,log=log))
          }
          )

pmcmc.internal <- function (object, Nmcmc = 1,
                            start = NULL,
                            pars = NULL, 
                            dprior.fun = NULL,
                            rw.sd = NULL,
                            Np = NULL,
                            hyperparams = NULL,
                            tol = 1e-17, max.fail = 0,
                            verbose = FALSE,
                            .ndone = 0, .prevcomp = NULL
                            ) {
  is.pmcmc <- is(object,"pmcmc")
  if (is.null(start)) {
    start <- coef(object)
    if (length(start)==0)
      stop("pmcmc error: ",sQuote("start")," must be specified if ",
           sQuote("coef(object)")," is NULL",
           call.=FALSE)
  } else if (length(start)==0)
    stop("pmcmc error: ",sQuote("start")," must be a named vector",call.=FALSE)
  start.names <- names(start)
  if (is.null(start.names))
    stop("pmcmc error: ",sQuote("start")," must be a named vector",call.=FALSE)

  if (is.null(rw.sd)) {
    if (is.pmcmc)
      rw.sd <- object@random.walk.sd
    else
      stop("pmcmc error: ",sQuote("rw.sd")," must be specified",call.=FALSE)
  }
  rw.names <- names(rw.sd)
  if (is.null(rw.names) || any(rw.sd<0))
    stop("pmcmc error: ",sQuote("rw.sd")," must be a named non-negative numerical vector",call.=FALSE)
  if (!all(rw.names%in%start.names))
    stop("pmcmc error: all the names of ",sQuote("rw.sd")," must be names of ",sQuote("start"),call.=FALSE)
  rw.names <- names(rw.sd[rw.sd>0])
  if (length(rw.names) == 0)
    stop("pmcmc error: ",sQuote("rw.sd")," must have one positive entry for each parameter to be estimated",call.=FALSE)

  if (is.null(pars)) {
    if (is.pmcmc) 
      pars <- object@pars
    else
      stop("pmcmc error: ",sQuote("pars")," must be specified",call.=FALSE)
  }
  if (length(pars)==0)
    stop("pmcmc error: at least one parameter must be estimated",call.=FALSE)
  if (
      !is.character(pars) ||
      !all(pars%in%start.names) ||
      !all(pars%in%rw.names)
      )
    stop(
         "pmcmc error: ",
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
            "pmcmc warning: the variable(s) ",
            paste(extra.rws,collapse=", "),
            " have positive random-walk SDs specified, but are not included in ",
            sQuote("pars"),
            ". These random walk SDs are ignored.",
            call.=FALSE
            )
  }
  rw.sd <- rw.sd[pars]
  rw.names <- names(rw.sd)

  if (is.null(dprior.fun)) {
    if (is.pmcmc) 
      dprior.fun <- object@dprior
    else
      stop("pmcmc error: ",sQuote("dprior")," must be specified",call.=FALSE)
  }

  if (is.null(Np)) {
    if (is.pmcmc) 
      Np <- object@Np
    else
      stop("pmcmc error: ",sQuote("Np")," must be specified",call.=FALSE)
  }
  Np <- as.integer(Np)
  if ((length(Np)!=1)||(Np < 1))
    stop("pmcmc error: ",sQuote("Np")," must be a positive integer",call.=FALSE)

  if (is.null(hyperparams)) {
    if (is.pmcmc) 
      hyperparams <- object@hyperparams
    else
      stop("pmcmc error: ",sQuote("hyperparams")," must be specified",call.=FALSE)
  }

  Nmcmc <- as.integer(Nmcmc)
  if (Nmcmc<0)
    stop("pmcmc error: ",sQuote("Nmcmc")," must be a positive integer",call.=FALSE)

  if (verbose) {
    cat("performing",Nmcmc,"PMCMC iteration(s) to estimate parameter(s)",
        paste(pars,collapse=", "))
    cat(" using random-walk with SD\n")
    print(rw.sd)
    cat("using",Np,"particles\n")
  }

  theta <- start

  conv.rec <- matrix(
                     data=NA,
                     nrow=Nmcmc+1,
                     ncol=length(theta)+3,
                     dimnames=list(
                       rownames=seq(from=0,to=Nmcmc,by=1),
                       colnames=c('loglik','log.prior','nfail',names(theta))
                       )
                     )

  if (!all(is.finite(theta[pars]))) {
    stop(
         sQuote("pmcmc"),
         " error: cannot estimate non-finite parameters: ",
         paste(
               pars[!is.finite(theta[pars])],
               collapse=","
               ),
         call.=FALSE
         )
  }

  obj <- new(
             "pmcmc",
             as(object,"pomp"),
             pars=pars,
             Nmcmc=0L,
             dprior=dprior.fun,
             Np=Np,
             hyperparams=hyperparams,
             random.walk.sd=rw.sd
             )

  if (.ndone==0) { ## compute prior and likelihood on initial parameter vector
    x <- try(
             pfilter.internal(
                              object=obj,
                              params=theta,
                              Np=Np,
                              tol=tol,
                              max.fail=max.fail,
                              pred.mean=FALSE,
                              pred.var=FALSE,
                              filter.mean=TRUE,
                              save.states=FALSE,
                              verbose=verbose
                              ),
             silent=FALSE
             )
    if (inherits(x,'try-error'))
      stop("pmcmc error: error in ",sQuote("pfilter"),call.=FALSE)
    x$log.prior <- dprior(object=obj,params=theta,log=TRUE)
  } else { ## has been computed previously
    x <- .prevcomp
  }
  conv.rec[1,names(theta)] <- theta
  conv.rec[1,c(1,2,3)] <- c(x$loglik,x$log.prior,x$nfail)

  for (n in seq_len(Nmcmc)) { # main loop

    theta.prop <- theta
    theta.prop[pars] <- rnorm(n=length(pars),mean=theta.prop[pars],sd=rw.sd)

    ## run the particle filter on the proposed new parameter values
    x.prop <- try(
                  pfilter.internal(
                                   object=obj,
                                   params=theta.prop,
                                   Np=Np,
                                   tol=tol,
                                   max.fail=max.fail,
                                   pred.mean=FALSE,
                                   pred.var=FALSE,
                                   filter.mean=TRUE,
                                   save.states=FALSE,
                                   verbose=verbose
                                   ),
                  silent=FALSE
                  )
    if (inherits(x.prop,'try-error'))
      stop("pmcmc error: error in ",sQuote("pfilter"),call.=FALSE)
    x.prop$log.prior <- dprior(object=obj,params=theta.prop,log=TRUE)

    ## PMCMC update rule (OK because proposal is symmetric)
    if (runif(1) < exp(x.prop$loglik-x$loglik+x.prop$log.prior-x$log.prior)) {
      theta <- theta.prop
      x <- x.prop
    }

    ## store a record of this iteration
    conv.rec[n+1,names(theta)] <- theta
    conv.rec[n+1,c(1,2,3)] <- c(x$loglik,x$log.prior,x$nfail)

    if (verbose) cat("PMCMC iteration ",n," of ",Nmcmc," completed\n")

  }

  coef(obj) <- theta

  if (Nmcmc>0) {
    obj@Nmcmc <- as.integer(Nmcmc)
    obj@filter.mean <- x$filter.mean
    obj@conv.rec <- conv.rec
    obj@eff.sample.size <- x$eff.sample.size
    obj@cond.loglik <- x$cond.loglik
    obj@loglik <- x$loglik
    obj@log.prior <- x$log.prior
  }

  obj
}

setMethod(
          "pmcmc",
          "pomp",
          function (
                    object, Nmcmc = 1,
                    start, pars, rw.sd,
                    dprior, Np, hyperparams,
                    tol = 1e-17, max.fail = 0, verbose = getOption("verbose")
                    )
          {
            if (missing(start)) start <- NULL
            if (missing(rw.sd))
              stop("pmcmc error: ",sQuote("rw.sd")," must be specified",call.=FALSE)
            if (missing(pars)) {
              pars <- names(rw.sd)[rw.sd>0]
            }
            if (missing(Np))
              stop("pmcmc error: ",sQuote("Np")," must be specified",call.=FALSE)
            if (missing(hyperparams))
              stop("pmcmc error: ",sQuote("hyperparams")," must be specified",call.=FALSE)
            if (missing(dprior)) {         # use default flat improper prior
              dprior <- function (params, hyperparams, log) {
                if (log) 0 else 1
              }
            } else {
              dprior <- match.fun(dprior)
              if (!all(c('params','hyperparams','log')%in%names(formals(dprior))))
                stop(
                     "pmcmc error: ",
                     sQuote("dprior"),
                     " must be a function of prototype ",
                     sQuote("dprior(params,hyperparams,log)"),
                     call.=FALSE
                     )
            }
              
            pmcmc.internal(
                           object=object,
                           Nmcmc=Nmcmc,
                           start=start,
                           pars=pars,
                           dprior.fun=dprior,
                           rw.sd=rw.sd,
                           Np=Np,
                           hyperparams=hyperparams,
                           tol=tol,
                           max.fail=max.fail,
                           verbose=verbose,
                           .ndone=0
                           )
          }
          )
          

setMethod(
          "pmcmc",
          "pmcmc",
          function (object, Nmcmc, ...)
          {
            if (missing(Nmcmc)) Nmcmc <- object@Nmcmc
            pmcmc.internal(object,Nmcmc=Nmcmc,...)
          }
          )

setMethod(
          'continue',
          'pmcmc',
          function (object, Nmcmc = 1, ...) {
            ndone <- object@Nmcmc
            obj <- pmcmc.internal(
                                  object=object,
                                  Nmcmc=Nmcmc,
                                  .ndone=ndone,
                                  .prevcomp=list(
                                    log.prior=object@log.prior,
                                    loglik=object@loglik,
                                    nfail=object@conv.rec[ndone+1,"nfail"],
                                    filter.mean=object@filter.mean,
                                    eff.sample.size=object@eff.sample.size,
                                    cond.loglik=object@cond.loglik
                                    ),
                                  ...
                                  )
            obj@conv.rec <- rbind(
                                  object@conv.rec[,colnames(obj@conv.rec)],
                                  obj@conv.rec[-1,]
                                  )
            obj@Nmcmc <- as.integer(ndone+Nmcmc)
            obj
          }
          )
