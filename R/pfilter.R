##' Particle filter
##'
##' A plain vanilla sequential Monte Carlo (particle filter) algorithm.
##' Resampling is performed at each observation.
##'
##' @name pfilter
##' @rdname pfilter
##' @aliases pfilter,ANY-method pfilter,missing-method
##' @author Aaron A. King
##' @family elementary algorithms
##' @family particle filter methods
##' @family full-information methods
##' @include pomp_class.R pomp.R rprocess_spec.R dmeasure_spec.R
##' @importFrom stats setNames
##' @inheritSection pomp Note for Windows users
##' @inheritParams pomp
##' @param Np the number of particles to use.
##' This may be specified as a single positive integer, in which case the same number of particles will be used at each timestep.
##' Alternatively, if one wishes the number of particles to vary across timesteps, one may specify \code{Np} either as a vector of positive integers of length \preformatted{length(time(object,t0=TRUE))} or as a function taking a positive integer argument.
##' In the latter case, \code{Np(k)} must be a single positive integer, representing the number of particles to be used at the \code{k}-th timestep:
##' \code{Np(0)} is the number of particles to use going from \code{timezero(object)} to \code{time(object)[1]},
##' \code{Np(1)}, from \code{timezero(object)} to \code{time(object)[1]},
##' and so on,
##' while when \code{T=length(time(object))}, \code{Np(T)} is the number of particles to sample at the end of the time-series.
##'
##' @param pred.mean logical; if \code{TRUE}, the prediction means are calculated for the state variables and parameters.
##'
##' @param pred.var logical; if \code{TRUE}, the prediction variances are calculated for the state variables and parameters.
##'
##' @param filter.mean logical; if \code{TRUE}, the filtering means are calculated for the state variables and parameters.
##'
##' @param filter.traj logical; if \code{TRUE}, a filtered trajectory is returned for the state variables and parameters.
##' See \code{\link{filter_traj}} for more information.
##'
##' @param save.states character;
##' If \code{save.states="unweighted"}, the state-vector for each unweighted particle at each time is saved.
##' If \code{save.states="weighted"}, the state-vector for each weighted particle at each time is saved, along with the corresponding weight.
##' If \code{save.states="no"}, information on the latent states is not saved.
##' \code{"FALSE"} is a synonym for \code{"no"} and \code{"TRUE"} is a synonym for \code{"unweighted"}.
##' To retrieve the saved states, applying \code{\link{saved.states}} to the result of the \code{pfilter} computation.
##'
##' @return
##' An object of class \sQuote{pfilterd_pomp}, which extends class \sQuote{pomp}.
##' Information can be extracted from this object using the methods documented below.
##' 
##' @section Methods:
##' \describe{
##' \item{\code{\link{logLik}}}{ the estimated log likelihood  }
##' \item{\code{\link{cond_logLik}}}{ the estimated conditional log likelihood }
##' \item{\code{\link{eff_sample_size}}}{
##' the (time-dependent) estimated effective sample size }
##' \item{\code{\link{pred_mean}}, \code{\link{pred_var}}}{ the mean and variance of the approximate prediction distribution }
##' \item{\code{\link{filter_mean}}}{ the mean of the filtering distribution }
##' \item{\code{\link{filter_traj}}}{
##'   retrieve one particle trajectory.
##'   Useful for building up the smoothing distribution.
##' }
##' \item{\code{\link{saved_states}}}{retrieve saved states}
##' \item{\code{\link{as.data.frame}}}{coerce to a data frame}
##' \item{\code{\link{plot}}}{diagnostic plots}
##' }
##'
##' @references
##'
##' \Arulampalam2002
##'
##' \Bhadra2016
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
    cond.logLik="numeric",
    saved.states="list",
    Np="integer",
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
    cond.logLik=numeric(0),
    saved.states=list(),
    Np=as.integer(NA),
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

##' @rdname pfilter
##' @export
setMethod(
  "pfilter",
  signature=signature(data="data.frame"),
  definition=function (
    data,
    Np, 
    params, rinit, rprocess, dmeasure,
    pred.mean = FALSE,
    pred.var = FALSE,
    filter.mean = FALSE,
    filter.traj = FALSE,
    save.states = c("no", "weighted", "unweighted", "FALSE", "TRUE"),
    ...,
    verbose = getOption("verbose", FALSE)) {

    tryCatch(
      pfilter_internal(
        data,
        Np=Np,
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

##' @rdname pfilter
##' @export
setMethod(
  "pfilter",
  signature=signature(data="pomp"),
  definition=function (
    data,
    Np,
    pred.mean = FALSE,
    pred.var = FALSE,
    filter.mean = FALSE,
    filter.traj = FALSE,
    save.states = c("no", "weighted", "unweighted", "FALSE", "TRUE"),
    ...,
    verbose = getOption("verbose", FALSE)) {

    tryCatch(
      pfilter_internal(
        data,
        Np=Np,
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

##' @rdname pfilter
##' @export
setMethod(
  "pfilter",
  signature=signature(data="pfilterd_pomp"),
  function (data, Np,
    ..., verbose = getOption("verbose", FALSE)) {

    if (missing(Np)) Np <- data@Np

    pfilter(as(data,"pomp"),Np=Np,...,verbose=verbose)

  }
)

pfilter_internal <- function (object, Np,
  pred.mean = FALSE, pred.var = FALSE, filter.mean = FALSE,
  filter.traj = FALSE, cooling, cooling.m,
  save.states = c("no", "weighted", "unweighted", "FALSE", "TRUE"),
  ..., .gnsi = TRUE, verbose = FALSE) {

  verbose <- as.logical(verbose)

  object <- pomp(object,...,verbose=verbose)

  if (undefined(object@rprocess) || undefined(object@dmeasure))
    pStop_(paste(sQuote(c("rprocess","dmeasure")),collapse=", ")," are needed basic components.")

  gnsi <- as.logical(.gnsi)
  pred.mean <- as.logical(pred.mean)
  pred.var <- as.logical(pred.var)
  filter.mean <- as.logical(filter.mean)
  filter.traj <- as.logical(filter.traj)
  save.states <- as.character(save.states)
  save.states <- match.arg(save.states)

  params <- coef(object)
  times <- time(object,t0=TRUE)
  ntimes <- length(times)-1

  Np <- np_check(Np,ntimes)

  pompLoad(object,verbose=verbose)
  on.exit(pompUnload(object,verbose=verbose))

  init.x <- rinit(object,params=params,nsim=Np[1L],.gnsi=gnsi)
  statenames <- rownames(init.x)
  nvars <- nrow(init.x)
  x <- init.x

  ## set up storage for saving samples from filtering distributions
  stsav <- save.states %in% c("unweighted","TRUE")
  wtsav <- save.states == "weighted"
  if (stsav || wtsav || filter.traj) {
    xparticles <- vector(mode="list",length=ntimes)
    if (wtsav) xweights <- xparticles
  }
  if (filter.traj) {
    pedigree <- vector(mode="list",length=ntimes+1)
  }

  loglik <- rep(NA,ntimes)
  eff.sample.size <- numeric(ntimes)

  ## set up storage for prediction means, variances, etc.
  if (pred.mean) {
    pred.m <- array(data=numeric(1),dim=c(nvars,ntimes),
      dimnames=list(name=statenames,time=NULL))
  } else {
    pred.m <- array(data=numeric(0),dim=c(0,0))
  }

  if (pred.var) {
    pred.v <- array(data=numeric(1),dim=c(nvars,ntimes),
      dimnames=list(name=statenames,time=NULL))
  } else {
    pred.v <- array(data=numeric(0),dim=c(0,0))
  }

  if (filter.mean) {
    filt.m <- array(data=numeric(1),dim=c(nvars,ntimes),
      dimnames=list(name=statenames,time=NULL))
  } else {
    filt.m <- array(data=numeric(0),dim=c(0,0))
  }

  if (filter.traj) {
    filt.t <- array(data=numeric(1),dim=c(nvars,1,ntimes+1),
      dimnames=list(name=statenames,rep=1,time=NULL))
  } else {
    filt.t <- array(data=numeric(0),dim=c(0,0,0))
  }

  for (nt in seq_len(ntimes)) { ## main loop

    ## advance the state variables according to the process model
    X <- rprocess(object,x0=x,t0=times[nt],times=times[nt+1],params=params,.gnsi=gnsi)

    if (pred.var) { ## check for nonfinite state variables and parameters
      problem.indices <- unique(which(!is.finite(X),arr.ind=TRUE)[,1L])
      nprob <- length(problem.indices)
      if (nprob > 0)
        pStop_("non-finite state variable",ngettext(nprob,"","s"),": ",
          paste(rownames(X)[problem.indices],collapse=', '))
    }

    ## determine the weights
    weights <- dmeasure(object,y=object@data[,nt,drop=FALSE],x=X,
      times=times[nt+1],params=params,log=TRUE,.gnsi=gnsi)
    gnsi <- FALSE

    ## store unweighted particles and their weights
    if (wtsav) {
      xparticles[[nt]] <- array(X,dim=dim(X)[c(1L,2L)],
        dimnames=list(name=statenames,.id=NULL))
      xweights[[nt]] <- array(weights,dim=length(weights),
        dimnames=list(.id=NULL))
    }

    ## compute prediction mean, prediction variance, filtering mean,
    ## effective sample size, log-likelihood.
    ## also do resampling.
    xx <- .Call(P_pfilter_computations,x=X,params=params,Np=Np[nt+1],
      predmean=pred.mean,predvar=pred.var,filtmean=filter.mean,
      trackancestry=filter.traj,doparRS=FALSE,weights=weights,
      wave=FALSE)

    ## the following is triggered by the first illegal weight value
    if (is.integer(xx)) {
      illegal_dmeasure_error(
        time=times[nt+1],
        loglik=weights[xx],
        datvals=object@data[,nt],
        states=X[,xx,1L],
        params=params
      )
    }

    loglik[nt] <- xx$loglik
    eff.sample.size[nt] <- xx$ess

    x <- xx$states
    params <- xx$params[,1L]

    if (pred.mean) pred.m[,nt] <- xx$pm
    if (pred.var) pred.v[,nt] <- xx$pv
    if (filter.mean) filt.m[,nt] <- xx$fm
    if (filter.traj) pedigree[[nt]] <- xx$ancestry

    if (stsav || filter.traj) {
      xparticles[[nt]] <- x
      dimnames(xparticles[[nt]]) <- list(name=statenames,.id=NULL)
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

  if (stsav) {
    stsav <- xparticles
  } else if (wtsav) {
    stsav <- list(states=xparticles,weights=xweights)
  } else {
    stsav <- list()
  }

  new(
    "pfilterd_pomp",
    object,
    pred.mean=pred.m,
    pred.var=pred.v,
    filter.mean=filt.m,
    filter.traj=filt.t,
    paramMatrix=array(data=numeric(0),dim=c(0,0)),
    eff.sample.size=eff.sample.size,
    cond.logLik=loglik,
    saved.states=stsav,
    Np=as.integer(Np),
    loglik=sum(loglik)
  )
}

illegal_dmeasure_error <- function (time, loglik, datvals, states, params) {
  showvals <- c(time=time,loglik=loglik,datvals,states,params)
  m1 <- formatC(names(showvals),preserve.width="common")
  m2 <- formatC(showvals,digits=6,width=12,format="g",preserve.width="common")
  pStop_(
    sQuote("dmeasure")," with log=TRUE returns illegal value.\n",
    "Log likelihood, data, states, and parameters are:\n",
    paste0(m1,": ",m2,collapse="\n")
  )
}

np_check <- function (Np, ntimes) {
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

  as.integer(Np)
}
