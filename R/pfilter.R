##' Particle filter
##'
##' A plain vanilla sequential Monte Carlo (particle filter) algorithm.
##' Resampling is performed at each observation.
##'
##' @name pfilter
##' @rdname pfilter
##' @aliases pfilter pfilter,ANY-method pfilter,missing-method
##' pfilterd_pomp-class pfilterd_pomp
##' @author Aaron A. King
##' @family elementary POMP methods
##' @family particle filter methods
##'
##' @include pomp_class.R pomp.R rprocess_spec.R dmeasure_spec.R
##' @importFrom stats setNames
##'
##' @inheritParams pomp
##'
##' @param Np the number of particles to use.
##' This may be specified as a single positive integer, in which case the same number of particles will be used at each timestep.
##' Alternatively, if one wishes the number of particles to vary across timesteps, one may specify \code{Np} either as a vector of positive integers of length \preformatted{length(time(object,t0=TRUE))} or as a function taking a positive integer argument.
##' In the latter case, \code{Np(k)} must be a single positive integer, representing the number of particles to be used at the \code{k}-th timestep:
##' \code{Np(0)} is the number of particles to use going from \code{timezero(object)} to \code{time(object)[1]},
##' \code{Np(1)}, from \code{timezero(object)} to \code{time(object)[1]},
##' and so on,
##' while when \code{T=length(time(object,t0=TRUE))}, \code{Np(T)} is the number of particles to sample at the end of the time-series.
##'
##' @param tol non-negative numeric scalar;
##' particles with likelihood less than \code{tol} are considered to be incompatible with the data.
##' See the section on \emph{Filtering Failures} for more information.
##' In a future release, this argument will be removed.
##'
##' @param max.fail integer; the maximum number of filtering failures allowed (see below).
##' If the number of filtering failures exceeds this number, execution will terminate with an error.
##' By default, \code{max.fail} is set to infinity, so no error can be triggered.
##' In a future release, this argument will be removed.
##'
##' @param pred.mean logical; if \code{TRUE}, the prediction means are calculated for the state variables and parameters.
##'
##' @param pred.var logical; if \code{TRUE}, the prediction variances are calculated for the state variables and parameters.
##'
##' @param filter.mean logical; if \code{TRUE}, the filtering means are calculated for the state variables and parameters.
##'
##' @param filter.traj logical; if \code{TRUE}, a filtered trajectory is returned for the state variables and parameters.
##' See \code{\link{filter.traj}} for more information.
##'
##' @param save.states logical.
##' If \code{save.states=TRUE}, the state-vector for each particle at each time is saved.
##'
##' @return
##' An object of class \sQuote{pfilterd_pomp}, which extends class \sQuote{pomp}.
##' @section Methods:
##' \describe{
##' \item{\code{\link{logLik}}}{ the estimated log likelihood  }
##' \item{\code{\link{cond.logLik}}}{ the estimated conditional log likelihood }
##' \item{\code{\link{eff.sample.size}}}{
##' the (time-dependent) estimated effective sample size }
##' \item{\code{\link{pred.mean}}, \code{\link{pred.var}}}{ the mean and variance of the approximate prediction distribution }
##' \item{\code{\link{filter.mean}}}{ the mean of the filtering distribution }
##' \item{\code{\link{filter.traj}}}{
##'   retrieve one particle trajectory.
##'   Useful for building up the smoothing distribution.
##' }
##' \item{\code{\link{as.data.frame}}}{ coerce to a data frame }
##' \item{\code{\link{plot}}}{diagnostic plots}
##' }
##'
##' @section Filtering failures:
##' If the degree of disagreement between model and data becomes sufficiently large, a \dQuote{filtering failure} results.
##' A filtering failure occurs when, at some time point, none of the \code{Np} particles is compatible with the data.
##' In particular, if the conditional likelihood of a particle at any time is below the tolerance value \code{tol}, then that particle is considered to be uninformative and its likelihood is taken to be zero.
##' A filtering failure occurs when this is the case for all particles.
##' A warning is generated when this occurs unless the cumulative number of failures exceeds \code{max.fail}, in which case an error is generated.
##'
##' @references
##'
##' \Arulampalam2002
##'
##' @example examples/pfilter.R
##'
NULL

setClass(
  "pfilterd_pomp",
  contains="pomp",
  slots=c(
    pred.mean="array",
    pred.var="array",
    filter.mean="array",
    filter.traj="array",
    paramMatrix="array",
    indices="vector",
    eff.sample.size="numeric",
    cond.loglik="numeric",
    saved.states="list",
    Np="integer",
    tol="numeric",
    nfail="integer",
    loglik="numeric"
  ),
  prototype=prototype(
    pred.mean=array(data=numeric(0),dim=c(0,0)),
    pred.var=array(data=numeric(0),dim=c(0,0)),
    filter.mean=array(data=numeric(0),dim=c(0,0)),
    filter.traj=array(data=numeric(0),dim=c(0,0,0)),
    paramMatrix=array(data=numeric(0),dim=c(0,0)),
    indices=integer(0),
    eff.sample.size=numeric(0),
    cond.loglik=numeric(0),
    saved.states=list(),
    Np=as.integer(NA),
    tol=as.double(NA),
    nfail=as.integer(NA),
    loglik=as.double(NA)
  )
)

setGeneric(
  "pfilter",
  function (data, ...)
    standardGeneric("pfilter")
)

setMethod(
  "pfilter",
  signature=signature(data="missing"),
  definition=function (...) {
    reqd_arg("pfilter","data")
  }
)

setMethod(
  "pfilter",
  signature=signature(data="ANY"),
  definition=function (data, ...) {
    undef_method("pfilter",data)
  }
)

##' @name pfilter-data.frame
##' @aliases pfilter,data.frame-method
##' @rdname pfilter
##' @export
setMethod(
  "pfilter",
  signature=signature(data="data.frame"),
  definition=function (
    data,
    Np, tol = 1e-17, max.fail = Inf,
    params, rinit, rprocess, dmeasure,
    pred.mean = FALSE,
    pred.var = FALSE,
    filter.mean = FALSE,
    filter.traj = FALSE,
    save.states = FALSE,
    ...,
    verbose = getOption("verbose", FALSE)) {

    tryCatch(
      pfilter.internal(
        data,
        Np=Np,
        tol=tol,
        max.fail=max.fail,
        pred.mean=pred.mean,
        pred.var=pred.var,
        filter.mean=filter.mean,
        filter.traj=filter.traj,
        save.states=save.states,
        rinit=rinit,
        rprocess=rprocess,
        dmeasure=dmeasure,
        params=params,
        ...,
        verbose=verbose
      ),
      error = function (e) pStop("pfilter",conditionMessage(e))
    )

  }
)

##' @name pfilter-pomp
##' @aliases pfilter,pomp-method
##' @rdname pfilter
##' @export
setMethod(
  "pfilter",
  signature=signature(data="pomp"),
  definition=function (
    data,
    Np, tol = 1e-17, max.fail = Inf,
    pred.mean = FALSE,
    pred.var = FALSE,
    filter.mean = FALSE,
    filter.traj = FALSE,
    save.states = FALSE,
    ...,
    verbose = getOption("verbose", FALSE)) {

    tryCatch(
      pfilter.internal(
        data,
        Np=Np,
        tol=tol,
        max.fail=max.fail,
        pred.mean=pred.mean,
        pred.var=pred.var,
        filter.mean=filter.mean,
        filter.traj=filter.traj,
        save.states=save.states,
        ...,
        verbose=verbose
      ),
      error = function (e) pStop("pfilter",conditionMessage(e))
    )

  }
)

##' @name pfilter-pfilterd_pomp
##' @aliases pfilter,pfilterd_pomp-method
##' @rdname pfilter
##' @export
setMethod(
  "pfilter",
  signature=signature(data="pfilterd_pomp"),
  function (data, Np, tol,
    ..., verbose = getOption("verbose", FALSE)) {

    if (missing(Np)) Np <- data@Np
    if (missing(tol)) tol <- data@tol

    pfilter(as(data,"pomp"),Np=Np,tol=tol,
      ...,verbose=verbose)

  }
)

pfilter.internal <- function (object, Np, tol, max.fail,
  pred.mean = FALSE, pred.var = FALSE, filter.mean = FALSE,
  filter.traj = FALSE, cooling, cooling.m, save.states = FALSE, ...,
  .gnsi = TRUE, verbose = FALSE) {

  verbose <- as.logical(verbose)

  object <- pomp(object,...,verbose=verbose)

  if (undefined(object@rprocess) || undefined(object@dmeasure))
    pStop_(paste(sQuote(c("rprocess","dmeasure")),collapse=", ")," are needed basic components.")

  tol <- as.numeric(tol)
  gnsi <- as.logical(.gnsi)
  pred.mean <- as.logical(pred.mean)
  pred.var <- as.logical(pred.var)
  filter.mean <- as.logical(filter.mean)
  filter.traj <- as.logical(filter.traj)
  save.states <- as.logical(save.states)

  params <- coef(object)
  times <- time(object,t0=TRUE)
  ntimes <- length(times)-1

  if (missing(Np) || is.null(Np)) {
    pStop_(sQuote("Np")," must be specified.")
  } else if (is.function(Np)) {
    Np <- tryCatch(
      vapply(seq.int(from=0,to=ntimes,by=1),Np,numeric(1)),
      error = function (e) {
        pStop_("if ",sQuote("Np")," is a function, it must return ",
          "a single positive integer.")
      }
    )
  } else if (!is.numeric(Np)) {
    pStop_(sQuote("Np")," must be a number, a vector of numbers, or a function.")
  }

  if (length(Np) == 1)
    Np <- rep(Np,times=ntimes+1)
  else if (length(Np) != (ntimes+1))
    pStop_(sQuote("Np")," must have length 1 or length ",ntimes+1,".")

  if (!all(is.finite(Np)) || any(Np <= 0))
    pStop_("number of particles, ",sQuote("Np"),", must be a positive integer.")

  Np <- as.integer(Np)

  if (length(tol) != 1 || !is.finite(tol) || tol < 0)
    pStop_(sQuote("tol")," should be a small nonnegative number.")

  if (tol != 0) {
    pWarn(
      "pfilter",
      "the ",sQuote("tol")," argument is deprecated and will be removed in a future release.\n",
      "Currently, the default value of ",sQuote("tol")," is 1e-17;\n",
      "in future releases, the value will be 0, and the option to choose otherwise will be removed."
    )
  }

  pompLoad(object,verbose=verbose)
  on.exit(pompUnload(object,verbose=verbose))

  init.x <- rinit(object,params=params,nsim=Np[1L],.gnsi=gnsi)
  statenames <- rownames(init.x)
  nvars <- nrow(init.x)
  x <- init.x

  ## set up storage for saving samples from filtering distributions
  if (save.states || filter.traj) {
    xparticles <- setNames(vector(mode="list",length=ntimes),time(object))
  }
  if (filter.traj) {
    pedigree <- vector(mode="list",length=ntimes+1)
  }

  loglik <- rep(NA,ntimes)
  eff.sample.size <- numeric(ntimes)
  nfail <- 0

  ## set up storage for prediction means, variances, etc.
  if (pred.mean) {
    pred.m <- array(data=numeric(1),dim=c(nvars,ntimes),
      dimnames=list(variable=statenames,time=time(object)))
  } else {
    pred.m <- array(data=numeric(0),dim=c(0,0))
  }

  if (pred.var) {
    pred.v <- array(data=numeric(1),dim=c(nvars,ntimes),
      dimnames=list(variable=statenames,time=time(object)))
  } else {
    pred.v <- array(data=numeric(0),dim=c(0,0))
  }

  if (filter.mean) {
    filt.m <- array(data=numeric(1),dim=c(nvars,ntimes),
      dimnames=list(variable=statenames,time=time(object)))
  } else {
    filt.m <- array(data=numeric(0),dim=c(0,0))
  }

  if (filter.traj) {
    filt.t <- array(data=numeric(1),dim=c(nvars,1,ntimes+1),
      dimnames=list(variable=statenames,rep=1,time=times))
  } else {
    filt.t <- array(data=numeric(0),dim=c(0,0,0))
  }

  for (nt in seq_len(ntimes)) { ## main loop

    ## advance the state variables according to the process model
    X <- rprocess(object,x0=x,t0=times[nt],times=times[nt+1],params=params,
      .gnsi=gnsi)

    if (pred.var) { ## check for nonfinite state variables and parameters
      problem.indices <- unique(which(!is.finite(X),arr.ind=TRUE)[,1L])
      nprob <- length(problem.indices)
      if (nprob > 0)
        pStop_("non-finite state variable",ngettext(nprob,"","s"),": ",
          paste(rownames(X)[problem.indices],collapse=', '))
    }

    ## determine the weights
    weights <- dmeasure(object,y=object@data[,nt,drop=FALSE],x=X,
      times=times[nt+1],params=params,log=FALSE,.gnsi=gnsi)

    gnsi <- FALSE

    ## compute prediction mean, prediction variance, filtering mean,
    ## effective sample size, log-likelihood
    ## also do resampling if filtering has not failed
    xx <- .Call(P_pfilter_computations,x=X,params=params,Np=Np[nt+1],
      predmean=pred.mean,predvar=pred.var,filtmean=filter.mean,
      trackancestry=filter.traj,doparRS=FALSE,weights=weights,
      wave=FALSE,tol=tol)

    ## the following is triggered by the first illegal weight value
    if (is.integer(xx)) {
      illegal_dmeasure_error(
        time=times[nt+1],
        lik=weights[xx],
        datvals=object@data[,nt],
        states=X[,xx,1L],
        params=params
      )
    }

    all.fail <- xx$fail
    loglik[nt] <- xx$loglik
    eff.sample.size[nt] <- xx$ess

    x <- xx$states
    params <- xx$params[,1L]

    if (pred.mean) pred.m[,nt] <- xx$pm
    if (pred.var) pred.v[,nt] <- xx$pv
    if (filter.mean) filt.m[,nt] <- xx$fm
    if (filter.traj) pedigree[[nt]] <- xx$ancestry

    if (all.fail) { ## all particles are lost
      nfail <- nfail+1
      if (verbose) message("filtering failure at time t = ",times[nt+1])
      if (nfail>max.fail) pStop_("too many filtering failures")
    }

    if (save.states || filter.traj) {
      xparticles[[nt]] <- x
      dimnames(xparticles[[nt]]) <- setNames(dimnames(xparticles[[nt]]),
        c("variable","rep"))
    }

    if (verbose && (nt%%5 == 0))
      cat("pfilter timestep",nt,"of",ntimes,"finished\n")

  } ## end of main loop

  if (filter.traj) { ## select a single trajectory
    b <- sample.int(n=length(weights),size=1L,replace=TRUE)
    filt.t[,1L,ntimes+1] <- xparticles[[ntimes]][,b]
    for (nt in seq.int(from=ntimes-1,to=1L,by=-1L)) {
      b <- pedigree[[nt+1]][b]
      filt.t[,1L,nt+1] <- xparticles[[nt]][,b]
    }
    if (times[2L] > times[1L]) {
      b <- pedigree[[1L]][b]
      filt.t[,1L,1L] <- init.x[,b]
    } else {
      filt.t <- filt.t[,,-1L,drop=FALSE]
    }
  }

  if (!save.states) xparticles <- list()

  if (nfail>0)
    pWarn("pfilter",nfail," filtering failure",ngettext(nfail,"","s")," occurred.")

  new(
    "pfilterd_pomp",
    object,
    pred.mean=pred.m,
    pred.var=pred.v,
    filter.mean=filt.m,
    filter.traj=filt.t,
    paramMatrix=array(data=numeric(0),dim=c(0,0)),
    eff.sample.size=eff.sample.size,
    cond.loglik=loglik,
    saved.states=xparticles,
    Np=as.integer(Np),
    tol=tol,
    nfail=as.integer(nfail),
    loglik=sum(loglik)
  )
}

illegal_dmeasure_error <- function (time, lik, datvals, states, params) {
  showvals <- c(time=time,lik=lik,datvals,states,params)
  m1 <- formatC(names(showvals),preserve.width="common")
  m2 <- formatC(showvals,digits=6,width=12,format="g",preserve.width="common")
  pStop_(
    sQuote("dmeasure")," returns illegal value.\n",
    "Likelihood, data, states, and parameters are:\n",
    paste0(m1,": ",m2,collapse="\n")
  )
}
