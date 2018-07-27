## define the abc class
setClass(
  'abc',
  contains='pomp',
  slots=c(
    pars = 'character',
    Nabc = 'integer',
    accepts = 'integer',
    probes='list',
    scale = 'numeric',
    epsilon = 'numeric',
    proposal = 'function',
    conv.rec = 'matrix'
  ),
  prototype=prototype(
    pars = character(0),
    Nabc = 0L,
    accepts = 0L,
    probes = list(),
    scale = numeric(0),
    epsilon = 1.0,
    proposal = function (...)
      stop("in ",sQuote("abc"),": proposal not specified",call.=FALSE),
    conv.rec=array(dim=c(0,0))
  )
)

setMethod(
  "abc",
  signature=signature(object="pomp"),
  definition=function (object, Nabc = 1,
    start, proposal,
    probes, scale, epsilon,
    verbose = getOption("verbose"),
    ...) {

    if (missing(start)) start <- coef(object)
    if (missing(proposal)) proposal <- NULL
    if (missing(probes)) probes <- NULL
    if (missing(scale)) scale <- NULL
    if (missing(epsilon)) epsilon <- NULL

    abc.internal(
      object=object,
      Nabc=Nabc,
      start=start,
      proposal=proposal,
      probes=probes,
      scale=scale,
      epsilon=epsilon,
      verbose=verbose,
      ...
    )
  }
)

setMethod(
  "abc",
  signature=signature(object="probed.pomp"),
  definition=function (object, probes,
    verbose = getOption("verbose"),
    ...) {

    if (missing(probes)) probes <- object@probes
    abc(object=as(object,"pomp"),probes=probes,...)

  }
)

setMethod(
  "abc",
  signature=signature(object="abc"),
  definition=function (object, Nabc,
    start, proposal,
    probes, scale, epsilon,
    verbose = getOption("verbose"),
    ...) {

    if (missing(Nabc)) Nabc <- object@Nabc
    if (missing(start)) start <- coef(object)
    if (missing(proposal)) proposal <- object@proposal
    if (missing(probes)) probes <- object@probes
    if (missing(scale)) scale <- object@scale
    if (missing(epsilon)) epsilon <- object@epsilon

    abc(
      object=as(object,"pomp"),
      Nabc=Nabc,
      start=start,
      proposal=proposal,
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
  definition=function (object, Nabc = 1, ...) {

    ndone <- object@Nabc
    accepts <- object@accepts

    obj <- abc(
      object=object,
      Nabc=Nabc,
      .ndone=ndone,
      .accepts=accepts,
      ...
    )

    obj@conv.rec <- rbind(
      object@conv.rec[,colnames(obj@conv.rec)],
      obj@conv.rec[-1,]
    )
    names(dimnames(obj@conv.rec)) <- c("iteration","variable")
    obj@Nabc <- as.integer(ndone+Nabc)
    obj@accepts <- as.integer(accepts+obj@accepts)

    obj
  }
)

abc.internal <- function (object, Nabc, start, proposal, probes, epsilon, scale,
  verbose = FALSE, .ndone = 0L, .accepts = 0L, .getnativesymbolinfo = TRUE,
  ...) {

  ep <- paste0("in ",sQuote("abc"),": ")

  object <- pomp(object,...)

  gnsi <- as.logical(.getnativesymbolinfo)
  .ndone <- as.integer(.ndone)
  .accepts <- as.integer(.accepts)
  scale <- as.numeric(scale)
  epsilon <- as.numeric(epsilon)
  epssq <- epsilon*epsilon
  verbose <- as.logical(verbose)

  Nabc <- as.integer(Nabc)
  if (!is.finite(Nabc) || Nabc < 0)
    stop(ep,sQuote("Nabc")," must be a positive integer.",call.=FALSE)

  if (is.list(start)) start <- unlist(start)
  start <- setNames(as.numeric(start),names(start))
  if (length(start)==0)
    stop(ep,"parameters must be specified.",call.=FALSE)
  start.names <- names(start)
  if (is.null(start.names) || !is.numeric(start))
    stop(ep,sQuote("start")," must be a named numeric vector.",call.=FALSE)

  if (is.null(proposal))
    stop(ep,sQuote("proposal")," must be specified",call.=FALSE)
  if (!is.function(proposal))
    stop(ep,sQuote("proposal")," must be a function.",call.=FALSE)

  if (is.null(probes))
    stop(ep,sQuote("probes")," must be specified",call.=FALSE)
  if (!is.list(probes)) probes <- list(probes)
  if (!all(sapply(probes,is.function)))
    stop(ep,sQuote("probes")," must be a function or a list of functions.",call.=FALSE)
  if (!all(sapply(probes,function(f)length(formals(f))==1)))
    stop(ep,"each probe must be a function of a single argument",call.=FALSE)

  if (length(scale)==0)
    stop(ep,sQuote("scale")," must be specified",call.=FALSE)

  if (length(epsilon)==0)
    stop(ep,"abc match criterion, ",sQuote("epsilon"),
      ", must be specified",call.=FALSE)

  pompLoad(object,verbose=verbose)

  ## test proposal distribution
  theta <- tryCatch(
    proposal(start,.n=0),
    error = function (e) {
      stop(ep,"in proposal function: ",conditionMessage(e),call.=FALSE)
    }
  )
  if (is.null(names(theta)) || !is.numeric(theta) || any(names(theta)==""))
    stop(ep,sQuote("proposal")," must return a named numeric vector.",call.=FALSE)


  if (verbose) {
    cat("performing",Nabc,"ABC iteration(s)\n")
  }

  theta <- start
  log.prior <- tryCatch(
    dprior(object,params=theta,log=TRUE,.getnativesymbolinfo=gnsi),
    error = function (e) {
      stop(ep,sQuote("dprior")," error: ",conditionMessage(e),call.=FALSE)
    }
  )
  if (!is.finite(log.prior))
    stop(ep,"inadmissible value of ",sQuote("dprior")," at start parameters.",
      call.=FALSE)
  ## we suppose that theta is a "match",
  ## which does the right thing for continue() and
  ## should have negligible effect unless doing many short calls to continue()

  conv.rec <- matrix(
    data=NA,
    nrow=Nabc+1,
    ncol=length(theta),
    dimnames=list(
      iteration=seq(from=0,to=Nabc,by=1),
      variable=names(theta)
    )
  )

  ## apply probes to data
  datval <- tryCatch(
    .Call(apply_probe_data,object,probes),
    error = function (e) {
      stop(ep,"in ",sQuote("apply_probe_data"),": ",conditionMessage(e),call.=FALSE)
    }
  )

  if (length(scale) != 1 && length(scale) != length(datval))
    stop(ep,sQuote("scale")," must have either length 1 or length equal to the",
      " number of probes (here, ",length(datval),").",call.=FALSE)

  conv.rec[1,names(theta)] <- theta

  for (n in seq_len(Nabc)) { # main loop

    theta.prop <- tryCatch(
      proposal(theta,.n=n+.ndone,.accepts=.accepts,verbose=verbose),
      error = function (e) {
        stop(ep,"in proposal function: ",conditionMessage(e),call.=FALSE)
      }
    )
    log.prior.prop <- tryCatch(
      dprior(object,params=theta.prop,log=TRUE,.getnativesymbolinfo=gnsi),
      error = function (e) {
        stop(ep,sQuote("dprior")," error: ",conditionMessage(e),call.=FALSE)
      }
    )

    if (is.finite(log.prior.prop) && runif(1) < exp(log.prior.prop-log.prior)) {

      ## compute the probes for the proposed new parameter values

      simval <- tryCatch(
        .Call(
          apply_probe_sim,
          object=object,
          nsim=1L,
          params=theta.prop,
          seed=NULL,
          probes=probes,
          datval=datval,
          gnsi=gnsi
        ),
        error = function (e) {
          stop(ep,"in ",sQuote("apply_probe_sim"),": ",conditionMessage(e),call.=FALSE)
        }
      )

      ## ABC update rule
      distance <- sum(((datval-simval)/scale)^2)
      if( (is.finite(distance)) && (distance<epssq) ){
        theta <- theta.prop
        log.prior <- log.prior.prop
        .accepts <- .accepts+1L
      }

      gnsi <- FALSE

    }

    ## store a record of this iteration
    conv.rec[n+1,names(theta)] <- theta
    if (verbose && (n%%5==0))
      cat("ABC iteration",n+.ndone,"of",Nabc+.ndone,
        "completed\nacceptance ratio:",
        round(.accepts/(n+.ndone),3),"\n")
  }

  pars <- apply(conv.rec,2,function(x)diff(range(x))>0)
  pars <- names(pars[pars])

  pompUnload(object,verbose=verbose)

  new(
    'abc',
    object,
    params=theta,
    pars=pars,
    Nabc=Nabc,
    accepts=.accepts,
    probes=probes,
    scale=scale,
    epsilon=epsilon,
    proposal=proposal,
    conv.rec=conv.rec
  )

}
