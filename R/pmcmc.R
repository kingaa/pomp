##' The particle Markov chain Metropolis-Hastings algorithm
##'
##' The Particle MCMC algorithm for estimating the parameters of a
##' partially-observed Markov process.  Running \code{pmcmc} causes a particle
##' random-walk Metropolis-Hastings Markov chain algorithm to run for the
##' specified number of proposals.
##'
##' @name pmcmc
##' @rdname pmcmc
##' @include pfilter.R proposals.R load.R continue.R
##' @aliases pmcmc,ANY-method pmcmc,missing-method
##' @author Edward L. Ionides, Aaron A. King, Sebastian Funk
##' @family estimation methods
##' @family particle filter methods
##' @family full-information methods
##' @family MCMC methods
##' @family Bayesian methods
##' 
##' @importFrom stats runif
##' @inheritParams pomp
##' @inheritParams pfilter
##' @param Nmcmc The number of PMCMC iterations to perform.
##' @param proposal optional function that draws from the proposal
##' distribution.  Currently, the proposal distribution must be symmetric for
##' proper inference: it is the user's responsibility to ensure that it is.
##' Several functions that construct appropriate proposal function are
##' provided: see \link[=proposals]{MCMC proposals} for more information.
##'
##' @section Methods:
##' The following can be applied to the output of a \code{pmcmc} operation:
##' \describe{
##' \item{\code{pmcmc}}{repeats the calculation, beginning with the last state}
##' \item{\code{\link{continue}}}{continues the \code{pmcmc} calculation}
##' \item{\code{plot}}{produces a series of diagnostic plots}
##' \item{\code{\link{filter.traj}}}{extracts a random sample from the smoothing distribution}
##' \item{\code{\link{traces}}}{produces an \code{\link[coda]{mcmc}} object, to which the various \pkg{coda} convergence diagnostics can be applied}
##' }
##'
##' @return An object of class \sQuote{pmcmcd_pomp}.
##'
##' @inheritSection pomp Note for Windows users
##' 
##' @section Re-running PMCMC Iterations:
##' To re-run a sequence of PMCMC
##' iterations, one can use the \code{pmcmc} method on a \sQuote{pmcmc} object.
##' By default, the same parameters used for the original PMCMC run are re-used
##' (except for \code{verbose}, the default of which is shown above).  If one
##' does specify additional arguments, these will override the defaults.
##'
##' @references
##'
##' \Andrieu2010
##' 
NULL

setClass(
  "pmcmcd_pomp",
  contains="pfilterd_pomp",
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
      pStop("pmcmc","proposal not specified."),
    traces=array(dim=c(0,0)),
    log.prior=numeric(0)
  )
)

setGeneric(
  "pmcmc",
  function (data, ...)
    standardGeneric("pmcmc")
)

setMethod(
  "pmcmc",
  signature=signature(data="missing"),
  definition=function (...) {
    reqd_arg("pmcmc","data")
  }
)

setMethod(
  "pmcmc",
  signature=signature(data="ANY"),
  definition=function (data, ...) {
    undef_method("pmcmc",data)
  }
)

##' @rdname pmcmc
##' @export
setMethod(
  "pmcmc",
  signature=signature(data="data.frame"),
  function (data,
    Nmcmc = 1, proposal,
    Np,
    params, rinit, rprocess, dmeasure, dprior,
    ..., verbose = getOption("verbose", FALSE)) {

    tryCatch(
      pmcmc.internal(
        data,
        Nmcmc=Nmcmc,
        proposal=proposal,
        Np=Np,
        params=params,
        rinit=rinit,
        rprocess=rprocess,
        dmeasure=dmeasure,
        dprior=dprior,
        ...,
        verbose=verbose
      ),
      error = function (e) pStop("pmcmc",conditionMessage(e))
    )

  }
)

##' @rdname pmcmc
##' @export
setMethod(
  "pmcmc",
  signature=signature(data="pomp"),
  function (data,
    Nmcmc = 1, proposal,
    Np, ...,
    verbose = getOption("verbose", FALSE)) {

    tryCatch(
      pmcmc.internal(
        data,
        Nmcmc=Nmcmc,
        proposal=proposal,
        Np=Np,
        ...,
        verbose=verbose
      ),
      error = function (e) pStop("pmcmc",conditionMessage(e))
    )

  }
)

##' @rdname pmcmc
##' @export
setMethod(
  "pmcmc",
  signature=signature(data="pfilterd_pomp"),
  function (data,
    Nmcmc = 1, proposal,
    Np, ...,
    verbose = getOption("verbose", FALSE)) {

    if (missing(Np)) Np <- data@Np
    
    pmcmc(
      as(data,"pomp"),
      Nmcmc=Nmcmc,
      proposal=proposal,
      Np=Np,
      ...,
      verbose=verbose
    )

  }
)

##' @rdname pmcmc
##' @export
setMethod(
  "pmcmc",
  signature=signature(data="pmcmcd_pomp"),
  function (data,
    Nmcmc, proposal,
    ..., verbose = getOption("verbose", FALSE)) {

    if (missing(Nmcmc)) Nmcmc <- data@Nmcmc
    if (missing(proposal)) proposal <- data@proposal

    pmcmc(
      as(data,"pfilterd_pomp"),
      Nmcmc=Nmcmc,
      proposal=proposal,
      ...,
      verbose=verbose
    )

  }
)

##' @rdname continue
##' @param Nmcmc positive integer; number of additional PMCMC iterations to perform
##' @export
setMethod(
  "continue",
  signature=signature(object="pmcmcd_pomp"),
  function (object, Nmcmc = 1, ...) {

    ndone <- object@Nmcmc
    accepts <- object@accepts

    obj <- pmcmc(
      object,
      Nmcmc=Nmcmc,
      ...,
      .ndone=ndone,
      .accepts=accepts,
      .prev.pfp=as(object,"pfilterd_pomp"),
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

pmcmc.internal <- function (object, Nmcmc, proposal, Np, ...,
  verbose, .ndone = 0L, .accepts = 0L, .prev.pfp = NULL, .prev.log.prior = NULL,
  .gnsi = TRUE) {

  verbose <- as.logical(verbose)

  object <- pomp(object,...,verbose=verbose)

  gnsi <- as.logical(.gnsi)
  .ndone <- as.integer(.ndone)
  .accepts <- as.integer(.accepts)

  if (missing(proposal) || is.null(proposal))
    pStop_(sQuote("proposal")," must be specified")
  proposal <- tryCatch(
    match.fun(proposal),
    error = function (e) {
      pStop_(sQuote("proposal")," must be a function: ",conditionMessage(e))
    }
  )

  pompLoad(object,verbose=verbose)
  on.exit(pompUnload(object,verbose=verbose))

  start <- coef(object)

  ## test proposal distribution
  theta <- tryCatch(
    proposal(start,.n=0L),
    error = function (e) {
      pStop_("error in proposal function: ",conditionMessage(e))
    }
  )
  if (is.null(names(theta)) || !is.numeric(theta) || any(names(theta)==""))
    pStop_(sQuote("proposal")," must return a named numeric vector.")

  theta <- start

  Nmcmc <- as.integer(Nmcmc)
  if (length(Nmcmc)!=1L || !is.finite(Nmcmc) || Nmcmc < 0L)
    pStop_(sQuote("Nmcmc")," must be a positive integer")

  if (verbose)
    cat("performing",Nmcmc,"PMCMC iteration(s) using",Np[1L],"particles\n")

  traces <- array(data=NA_real_,dim=c(Nmcmc+1L,length(theta)+2L),
    dimnames=list(
      iteration=seq(from=0L,to=Nmcmc,by=1L),
      variable=c("loglik","log.prior",names(theta))))

  if (.ndone==0L) { ## compute prior and likelihood on initial parameter vector

    log.prior <- dprior(object,params=theta,log=TRUE,.gnsi=gnsi)

    if (!is.finite(log.prior))
      pStop_("non-finite log prior at starting parameters.")

    pfp <- pfilter(
      object,
      params=theta,
      Np=Np,
      filter.traj=TRUE,
      .gnsi=gnsi,
      verbose=verbose
    )

    if (!is.finite(pfp@loglik))
      pStop_("non-finite log likelihood at starting parameters.")

  } else { ## has been computed previously

    pfp <- .prev.pfp
    log.prior <- .prev.log.prior
    pfp@filter.traj <- pfp@filter.traj[,.ndone,,drop=FALSE]

  }

  traces[1L,names(theta)] <- theta
  traces[1L,c(1L,2L)] <- c(pfp@loglik,log.prior)

  filt.t <- array(
    data=NA_real_,
    dim=replace(dim(pfp@filter.traj),2L,Nmcmc),
    dimnames=replace(
      dimnames(pfp@filter.traj),2L,list(as.character(seq_len(Nmcmc)))
    )
  )

  for (n in seq_len(Nmcmc)) { # main loop

    theta.prop <- tryCatch(
      proposal(theta,.n=n+.ndone,.accepts=.accepts,verbose=verbose),
      error = function (e) {
        pStop_("error in proposal function: ",conditionMessage(e))
      }
    )

    ## compute log prior
    log.prior.prop <- dprior(object,params=theta.prop,log=TRUE,.gnsi=gnsi)

    if (is.finite(log.prior.prop)) {

      ## run the particle filter on the proposed new parameter values
      pfp.prop <- pfilter(
        pfp,
        params=theta.prop,
        Np=Np,
        filter.traj=TRUE,
        verbose=verbose,
        .gnsi=gnsi
      )

      ## PMCMC update rule (OK because proposal is symmetric)
      alpha <- exp(pfp.prop@loglik+log.prior.prop-pfp@loglik-log.prior)

      if (!is.finite(alpha))
        pWarn_("non-finite log likelihood or log prior encountered.") #nocov

      if (is.finite(alpha) && runif(1) < alpha) {
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
    traces[n+1,c(1L,2L)] <- c(pfp@loglik,log.prior)

    if (verbose) cat("PMCMC iteration",n+.ndone,"of",Nmcmc+.ndone,
      "completed\nacceptance ratio:",
      round(.accepts/(n+.ndone),3),"\n")

    gnsi <- FALSE

  }

  pars <- apply(traces,2,function(x)diff(range(x))>0)
  pars <- setdiff(names(pars[pars]),c("loglik","log.prior"))

  new(
    "pmcmcd_pomp",
    pfp,
    params=theta,
    pars=pars,
    Nmcmc=Nmcmc,
    accepts=.accepts,
    proposal=proposal,
    traces=traces,
    log.prior=log.prior,
    filter.traj=filt.t
  )

}
