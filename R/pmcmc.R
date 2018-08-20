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
##' @aliases pmcmc pmcmc,ANY-method pmcmc,missing-method
##' @author Edward L. Ionides, Aaron A. King, Sebastian Funk
##' @family particle filter methods
##' @seealso \link[=proposals]{MCMC proposals}
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
##' @return An object of class \sQuote{pmcmcd_pomp}.
##'
##' @section Re-running PMCMC Iterations:
##' To re-run a sequence of PMCMC
##' iterations, one can use the \code{pmcmc} method on a \sQuote{pmcmc} object.
##' By default, the same parameters used for the original PMCMC run are re-used
##' (except for \code{tol}, \code{max.fail}, and \code{verbose}, the defaults
##' of which are shown above).  If one does specify additional arguments, these
##' will override the defaults.
##'
##' @references
##' C. Andrieu, A. Doucet, and R. Holenstein (2010)
##' Particle Markov chain Monte Carlo methods.
##' Journal of the Royal Statistical Society, Series B, 72: 269–342.
##'
##' C. Andrieu and G.O. Roberts (2009)
##' The pseudo-marginal approach for computation
##' Annals of Statistics, 37:697-725.
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
      stop("in ",sQuote("pmcmc"),": proposal not specified",call.=FALSE),
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
    stop("in ",sQuote("pmcmc"),": ",sQuote("data")," is a required argument",call.=FALSE)
  }
)

setMethod(
  "pmcmc",
  signature=signature(data="ANY"),
  definition=function (data, ...) {
    stop(sQuote("pmcmc")," is not defined when ",sQuote("data")," is of class ",sQuote(class(data)),call.=FALSE)
  }
)

##' @name pmcmc-data.frame
##' @aliases pmcmc,data.frame-method
##' @rdname pmcmc
##' @export
setMethod(
  "pmcmc",
  signature=signature(data="data.frame"),
  function (data, Nmcmc = 1, params, rinit, rprocess, dmeasure, dprior,
    proposal, Np, tol = 1e-17, max.fail = Inf, ...,
    verbose = getOption("verbose", FALSE)) {

    object <- pomp(data,params=params,rinit=rinit,rprocess=rprocess,
      dmeasure=dmeasure,dprior=dprior,...,verbose=verbose)

    pmcmc(
      object,
      Nmcmc=Nmcmc,
      proposal=proposal,
      Np=Np,
      tol=tol,
      max.fail=max.fail,
      verbose=verbose
    )

  }
)

##' @name pmcmc-pomp
##' @aliases pmcmc,pomp-method
##' @rdname pmcmc
##' @export
setMethod(
  "pmcmc",
  signature=signature(data="pomp"),
  function (data, Nmcmc = 1, proposal, Np, tol = 1e-17,
    max.fail = Inf, ..., verbose = getOption("verbose", FALSE)) {

    pmcmc.internal(
      data,
      Nmcmc=Nmcmc,
      proposal=proposal,
      Np=Np,
      tol=tol,
      max.fail=max.fail,
      verbose=verbose,
      ...
    )

  }
)

##' @name pmcmc-pfilterd_pomp
##' @aliases pmcmc,pfilterd_pomp-method
##' @rdname pmcmc
##' @export
setMethod(
  "pmcmc",
  signature=signature(data="pfilterd_pomp"),
  function (data, Nmcmc = 1, Np, tol, max.fail=Inf, ...,
    verbose = getOption("verbose", FALSE)) {

    if (missing(Np)) Np <- data@Np
    if (missing(tol)) tol <- data@tol

    pmcmc(as(data,"pomp"),Nmcmc=Nmcmc,Np=Np,tol=tol,...)
  }
)

##' @name pmcmc-pmcmcd_pomp
##' @aliases pmcmc,pmcmcd_pomp-method
##' @rdname pmcmc
##' @export
setMethod(
  "pmcmc",
  signature=signature(data="pmcmcd_pomp"),
  function (data, Nmcmc, proposal, Np, tol, max.fail = Inf,
    ..., verbose = getOption("verbose", FALSE)) {

    if (missing(Nmcmc)) Nmcmc <- data@Nmcmc
    if (missing(proposal)) proposal <- data@proposal
    if (missing(Np)) Np <- data@Np
    if (missing(tol)) tol <- data@tol

    pmcmc(as(data,"pomp"),Nmcmc=Nmcmc,proposal=proposal,
      Np=Np,tol=tol,max.fail=max.fail,verbose=verbose,...)
  }
)

##' @name continue-pmcmcd_pomp
##' @aliases continue,pmcmcd_pomp-method
##' @rdname continue
##'
##' @param Nmcmc positive integer; number of additional iterations to perform
##'
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

pmcmc.internal <- function (object, Nmcmc, proposal, Np, tol, max.fail, ...,
  verbose, .ndone = 0L, .accepts = 0L, .prev.pfp = NULL, .prev.log.prior = NULL,
  .getnativesymbolinfo = TRUE) {

  ep <- paste0("in ",sQuote("pmcmc"),": ")

  gnsi <- as.logical(.getnativesymbolinfo)
  ntimes <- length(time(object))
  verbose <- as.logical(verbose)
  .ndone <- as.integer(.ndone)
  .accepts <- as.integer(.accepts)

  object <- pomp(object,...,verbose=verbose)

  if (missing(Np) || is.null(Np)) {
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
    stop(ep,sQuote("Np")," must have length 1 or length ",
      sQuote("length(time(object))+1"),".",call.=FALSE)

  if (!all(is.finite(Np)) || any(Np <= 0))
    stop(ep,"number of particles, ",sQuote("Np"),
      ", must be a positive integer.",call.=FALSE)
  Np <- as.integer(Np)

  if (missing(proposal) || is.null(proposal))
    stop(ep,sQuote("proposal")," must be specified",call.=FALSE)
  if (!is.function(proposal))
    stop(ep,sQuote("proposal")," must be a function",call.=FALSE)

  pompLoad(object,verbose=verbose)
  on.exit(pompUnload(object,verbose=verbose))

  start <- coef(object)

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
        object,
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

  new(
    "pmcmcd_pomp",
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
