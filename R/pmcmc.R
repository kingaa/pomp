## define the pmcmc class
setClass(
  "pmcmc",
  contains="pfilterd.pomp",
  slots=c(
    pars = "character",
    Nmcmc = "integer",
    accepts = "integer",
    proposal = "function",
    traces = "array",
    log.prior = "numeric"
  ),
  prototype=prototype(
    pars = character(0),
    Nmcmc = 0L,
    accepts = 0L,
    proposal = function (...)
      stop("in ",sQuote("pmcmc"),": proposal not specified",call.=FALSE),
    traces=array(dim=c(0,0)),
    log.prior=numeric(0)
  )
)

setMethod(
  "pmcmc",
  signature=signature(object="pomp"),
  function (object, Nmcmc = 1,
    start, proposal, Np,
    tol = 1e-17, max.fail = Inf,
    verbose = getOption("verbose"),
    ...) {

    if (missing(start)) start <- coef(object)
    if (missing(proposal)) proposal <- NULL

    pmcmc.internal(
      object=object,
      Nmcmc=Nmcmc,
      start=start,
      proposal=proposal,
      Np=Np,
      tol=tol,
      max.fail=max.fail,
      verbose=verbose,
      ...
    )
  }
)

setMethod(
  "pmcmc",
  signature=signature(object="pfilterd.pomp"),
  function (object, Nmcmc = 1, Np, tol, ...) {

    if (missing(Np)) Np <- object@Np
    if (missing(tol)) tol <- object@tol

    pmcmc(as(object,"pomp"),Nmcmc=Nmcmc,Np=Np,tol=tol,...)
  }
)

setMethod(
  "pmcmc",
  signature=signature(object="pmcmc"),
  function (object, Nmcmc, start, proposal,
    Np, tol, max.fail = Inf,
    verbose = getOption("verbose"),
    ...) {

    if (missing(Nmcmc)) Nmcmc <- object@Nmcmc
    if (missing(start)) start <- coef(object)
    if (missing(proposal)) proposal <- object@proposal
    if (missing(Np)) Np <- object@Np
    if (missing(tol)) tol <- object@tol

    pmcmc(as(object,"pomp"),Nmcmc=Nmcmc,start=start,proposal=proposal,
      Np=Np,tol=tol,max.fail=max.fail,verbose=verbose,...)
  }
)

setMethod(
  "continue",
  signature=signature(object="pmcmc"),
  function (object, Nmcmc = 1, ...) {

    ndone <- object@Nmcmc
    accepts <- object@accepts

    obj <- pmcmc(
      object=object,
      Nmcmc=Nmcmc,
      ...,
      .ndone=ndone,
      .accepts=accepts,
      .prev.pfp=as(object,"pfilterd.pomp"),
      .prev.log.prior=object@log.prior
    )

    obj@traces <- rbind(
      object@traces[,colnames(obj@traces)],
      obj@traces[-1,]
    )
    names(dimnames(obj@traces)) <- c("iteration","variable")
    ft <- array(dim=replace(dim(obj@filter.traj),2L,ndone+Nmcmc),
      dimnames=replace(dimnames(obj@filter.traj),2L,
        list(seq_len(ndone+Nmcmc))))
    ft[,seq_len(ndone),] <- object@filter.traj
    ft[,ndone+seq_len(Nmcmc),] <- obj@filter.traj
    obj@filter.traj <- ft
    obj@Nmcmc <- as.integer(ndone+Nmcmc)
    obj@accepts <- as.integer(accepts+obj@accepts)

    obj
  }
)

pmcmc.internal <- function (object, Nmcmc,
  start, proposal, Np, tol, max.fail, ...,
  verbose,
  .ndone = 0L, .accepts = 0L,
  .prev.pfp = NULL, .prev.log.prior = NULL,
  .getnativesymbolinfo = TRUE) {

  object <- pomp(object,...)

  gnsi <- as.logical(.getnativesymbolinfo)
  verbose <- as.logical(verbose)
  .ndone <- as.integer(.ndone)
  .accepts <- as.integer(.accepts)
  ntimes <- length(time(object))

  ep <- paste0("in ",sQuote("pmcmc"),": ")

  if (is.list(start)) start <- unlist(start)
  if (is.null(start)) start <- numeric(0)
  if (length(start)==0)
    stop(ep,sQuote("start")," must be specified",call.=FALSE)
  if (!is.numeric(start) || is.null(names(start)))
    stop(ep,sQuote("start")," must be a named numeric vector",call.=FALSE)

  if (is.null(Np)) {
    stop(ep,sQuote("Np")," must be specified.",call.=FALSE)
  } else if (is.function(Np)) {
    Np <- tryCatch(
      vapply(seq.int(from=0,to=ntimes,by=1),Np,numeric(1)),
      error = function (e) {
        stop(ep,"if ",sQuote("Np")," is a function, ",
          "it must return a single positive integer.",call.=FALSE)
      }
    )
  } else if (!is.numeric(Np)) {
    stop(ep,sQuote("Np")," must be a number, a vector of numbers, ",
      "or a function.",call.=FALSE)
  }

  if (length(Np) == 1)
    Np <- rep(Np,times=ntimes+1)
  else if (length(Np) != (ntimes+1))
    stop(ep,sQuote("Np")," must have length 1 or length ",ntimes+1,".",call.=FALSE)

  if (!all(is.finite(Np)) || any(Np <= 0))
    stop(ep,"number of particles, ",sQuote("Np"),
      ", must be a positive integer.",call.=FALSE)
  Np <- as.integer(Np)

  if (is.null(proposal))
    stop(ep,sQuote("proposal")," must be specified",call.=FALSE)
  if (!is.function(proposal))
    stop(ep,sQuote("proposal")," must be a function",call.=FALSE)

  pompLoad(object,verbose=verbose)

  ## test proposal distribution
  theta <- tryCatch(
    proposal(start,.n=0),
    error = function (e) {
      stop(ep,"error in proposal function: ",conditionMessage(e),call.=FALSE)
    }
  )
  if (is.null(names(theta)) || !is.numeric(theta) || any(names(theta)==""))
    stop(ep,sQuote("proposal")," must return a named numeric vector",call.=FALSE)

  Nmcmc <- as.integer(Nmcmc)
  if (length(Nmcmc)!=1 || !is.finite(Nmcmc) || Nmcmc < 0)
    stop(ep,sQuote("Nmcmc")," must be a positive integer",call.=FALSE)

  if (verbose)
    cat("performing",Nmcmc,"PMCMC iteration(s) using",Np[1L],"particles\n")

  traces <- matrix(
    data=NA,
    nrow=Nmcmc+1,
    ncol=length(theta)+3,
    dimnames=list(
      iteration=seq(from=0,to=Nmcmc,by=1),
      variable=c("loglik","log.prior","nfail",names(theta))
    )
  )

  if (.ndone==0L) { ## compute prior and likelihood on initial parameter vector
    pfp <- tryCatch(
      pfilter.internal(
        object=object,
        params=theta,
        Np=Np,
        tol=tol,
        max.fail=max.fail,
        pred.mean=FALSE,
        pred.var=FALSE,
        filter.mean=TRUE,
        filter.traj=TRUE,
        save.states=FALSE,
        save.params=FALSE,
        verbose=FALSE,
        .getnativesymbolinfo=gnsi
      ),
      error = function (e) {
        stop(ep,conditionMessage(e),call.=FALSE) # nocov
      }
    )
    log.prior <- tryCatch(
      dprior(object,params=theta,log=TRUE,.getnativesymbolinfo=gnsi),
      error = function (e) {
        stop(ep,sQuote("dprior")," error: ",conditionMessage(e),call.=FALSE)
      }
    )
    gnsi <- FALSE
  } else { ## has been computed previously
    pfp <- .prev.pfp
    log.prior <- .prev.log.prior
    pfp@filter.traj <- pfp@filter.traj[,.ndone,,drop=FALSE]
  }
  traces[1,names(theta)] <- theta
  traces[1,c(1,2,3)] <- c(pfp@loglik,log.prior,pfp@nfail)

  filt.t <- array(
    data=0,
    dim=replace(dim(pfp@filter.traj),2L,Nmcmc),
    dimnames=replace(dimnames(pfp@filter.traj),2L,
      list(as.character(seq_len(Nmcmc))))
  )

  for (n in seq_len(Nmcmc)) { # main loop

    theta.prop <- tryCatch(
      proposal(theta,.n=n+.ndone,.accepts=.accepts,verbose=verbose),
      error = function (e) {
        stop(ep,"error in proposal function: ",conditionMessage(e),call.=FALSE)
      }
    )

    ## compute log prior
    log.prior.prop <- tryCatch(
      dprior(object,params=theta.prop,log=TRUE,.getnativesymbolinfo=gnsi),
      error = function (e) {
        stop(ep,sQuote("dprior")," error: ",conditionMessage(e),call.=FALSE)
      }
    )

    if (is.finite(log.prior.prop)) {

      ## run the particle filter on the proposed new parameter values
      pfp.prop <- tryCatch(
        pfilter.internal(
          object=pfp,
          params=theta.prop,
          Np=Np,
          tol=tol,
          max.fail=max.fail,
          pred.mean=FALSE,
          pred.var=FALSE,
          filter.mean=TRUE,
          filter.traj=TRUE,
          save.states=FALSE,
          save.params=FALSE,
          verbose=FALSE,
          .getnativesymbolinfo=gnsi
        ),
        error = function (e) {
          stop(ep,conditionMessage(e),call.=FALSE) # nocov
        }
      )
      gnsi <- FALSE

      ## PMCMC update rule (OK because proposal is symmetric)
      alpha <- exp(pfp.prop@loglik+log.prior.prop-pfp@loglik-log.prior)
      if (runif(1) < alpha) {
        pfp <- pfp.prop
        theta <- theta.prop
        log.prior <- log.prior.prop
        .accepts <- .accepts+1L
      }
    }

    ## add filtered trajectory to the store
    filt.t[,n,] <- pfp@filter.traj[,1L,]

    ## store a record of this iteration
    traces[n+1,names(theta)] <- theta
    traces[n+1,c(1,2,3)] <- c(pfp@loglik,log.prior,pfp@nfail)

    if (verbose) cat("PMCMC iteration",n+.ndone,"of",Nmcmc+.ndone,
      "completed\nacceptance ratio:",
      round(.accepts/(n+.ndone),3),"\n")

  }

  pars <- apply(traces,2,function(x)diff(range(x))>0)
  pars <- setdiff(names(pars[pars]),c("loglik","log.prior","nfail"))

  pompUnload(object,verbose=verbose)

  new(
    "pmcmc",
    pfp,
    params=theta,
    pars=pars,
    Nmcmc=Nmcmc,
    accepts=.accepts,
    proposal=proposal,
    Np=Np,
    tol=tol,
    traces=traces,
    log.prior=log.prior,
    filter.traj=filt.t
  )
}
