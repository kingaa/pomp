## define the abc class
setClass(
         'abc',
         contains='pomp',
         slots=c(
           pars = 'character',
           Nabc = 'integer',
           probes='list',
           scale = 'numeric',
           epsilon = 'numeric',
           random.walk.sd = 'numeric',
           conv.rec = 'matrix'
           )
         )

abc.internal <- function (object, Nabc,
                          start, pars,
                          rw.sd, probes,
                          epsilon, scale,
                          verbose,
                          .ndone = 0L,
                          .getnativesymbolinfo = TRUE,
                          ...) {

  object <- as(object,'pomp')
  gnsi <- as.logical(.getnativesymbolinfo)
  Nabc <- as.integer(Nabc)
  epsilon <- as.numeric(epsilon)
  epssq <- epsilon*epsilon

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
         " must be a subset of ",
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
  log.prior <- dprior(object,params=theta,log=TRUE,.getnativesymbolinfo=gnsi)
  if (!is.finite(log.prior))
    stop("inadmissible value of ",sQuote("dprior"),
         " at parameters ",sQuote("start"))
  ## we suppose that theta is a "match",
  ## which does the right thing for continue() and
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

  ## apply probes to data
  datval <- try(.Call(apply_probe_data,object,probes),silent=FALSE)
  if (inherits(datval,'try-error'))
    stop("abc error: error in ",sQuote("apply_probe_data"),call.=FALSE)

  conv.rec[1,names(theta)] <- theta

  for (n in seq_len(Nabc)) { # main loop

    theta.prop <- theta
    theta.prop[pars] <- rnorm(n=length(pars),mean=theta[pars],sd=rw.sd)
    log.prior.prop <- dprior(object,params=theta.prop,log=TRUE)

    if (is.finite(log.prior.prop) &&
        runif(1) < exp(log.prior.prop-log.prior)) {

      ## compute the probes for the proposed new parameter values

      simval <- try(
                    .Call(
                          apply_probe_sim,
                          object=object,
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
      if( (is.finite(distance)) && (distance<epssq) ){ 
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
      object,
      params=theta,
      pars=pars,
      Nabc=Nabc,
      probes=probes,
      scale=scale,
      epsilon=epsilon,
      random.walk.sd=rw.sd,
      conv.rec=conv.rec
      )

}

setMethod(
          "abc",
          signature=signature(object="pomp"),
          function (object, Nabc = 1,
                    start, pars, rw.sd,
                    probes, scale, epsilon,
                    verbose = getOption("verbose"),
                    ...) {
            
            if (missing(start))
              start <- coef(object)

            if (missing(rw.sd))
              stop("abc error: ",sQuote("rw.sd")," must be specified",
                   call.=FALSE)

            if (missing(pars))
              pars <- names(rw.sd)[rw.sd>0]

            if (missing(probes))
              stop("abc error: ",sQuote("probes")," must be specified",
                   call.=FALSE)

            if (missing(scale))
              stop("abc error: ",sQuote("scale")," must be specified",
                   call.=FALSE)

            if (missing(epsilon))
              stop("abc error: abc match criterion, ",sQuote("epsilon"),
                   ", must be specified",call.=FALSE)

            abc.internal(
                         object=object,
                         Nabc=Nabc,
                         start=start,
                         pars=pars,
                         probes=probes,
                         scale=scale,
                         epsilon=epsilon,
                         rw.sd=rw.sd,
                         verbose=verbose
                         )
          }
          )

setMethod(
          "abc",
          signature=signature(object="probed.pomp"),
          function (object, probes,
                    verbose = getOption("verbose"),
                    ...) {

            if (missing(probes)) probes <- object@probes
            f <- selectMethod("abc","pomp")
            f(
              object=object,
              probes=probes,
              ...
              )
          }
          )

setMethod(
          "abc",
          signature=signature(object="abc"),
          function (object, Nabc,
                    start, pars, rw.sd,
                    probes, scale, epsilon,
                    verbose = getOption("verbose"),
                    ...) {

            if (missing(Nabc)) Nabc <- object@Nabc
            if (missing(start)) start <- coef(object)
            if (missing(pars)) pars <- object@pars
            if (missing(rw.sd)) rw.sd <- object@random.walk.sd
            if (missing(probes)) probes <- object@probes
            if (missing(scale)) scale <- object@scale
            if (missing(epsilon)) epsilon <- object@epsilon

            f <- selectMethod("abc","pomp")

            f(
              object=object,
              Nabc=Nabc,
              start=start,
              pars=pars,
              rw.sd=rw.sd,
              probes=probes,
              scale=scale,
              epsilon=epsilon,
              verbose=verbose,
              ...
              )
          }
          )

setMethod(
          'continue',
          signature=signature(object='abc'),
          function (object, Nabc = 1, ...) {

            ndone <- object@Nabc
            f <- selectMethod("abc","abc")
            
            obj <- f(
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
