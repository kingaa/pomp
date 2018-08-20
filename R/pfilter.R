##' Particle filter
##'
##' A plain vanilla sequential Monte Carlo (particle filter) algorithm.
##' Resampling is performed at each observation.
##'
##' @name pfilter
##' @rdname pfilter
##' @aliases pfilter,ANY-method pfilter,missing-method pfilter
##' @author Aaron A. King
##' @family elementary POMP methods
##' @family particle filter methods
##'
##' @include pomp_class.R pomp.R rprocess_spec.R dmeasure_spec.R
##' @importFrom stats setNames
##'
##' @inheritParams pomp
##' @param Np the number of particles to use.
##' This may be specified as a single positive integer, in which case the same number of particles will be used at each timestep.
##' Alternatively, if one wishes the number of particles to vary across timesteps, one may specify \code{Np} either as a vector of positive integers of length \preformatted{length(time(object,t0=TRUE))} or as a function taking a positive integer argument.
##' In the latter case, \code{Np(k)} must be a single positive integer, representing the number of particles to be used at the \code{k}-th timestep:
##' \code{Np(0)} is the number of particles to use going from \code{timezero(object)} to \code{time(object)[1]},
##' \code{Np(1)}, from \code{timezero(object)} to \code{time(object)[1]},
##' and so on,
##' while when \code{T=length(time(object,t0=TRUE))}, \code{Np(T)} is the number of particles to sample at the end of the time-series.
##' When \code{object} is of class \sQuote{mif}, this is by default the same number of particles used in the \code{mif} iterations.
##'
##' One should omit \code{Np} if \code{params} is a matrix of parameters, with one column for each particle.  In this case, obviously, the number of particles is \code{ncol(params)}.
##' @param tol positive numeric scalar;
##' particles with likelihood less than \code{tol} are considered to be incompatible with the data.
##' See the section on \emph{Filtering Failures} below for more information.
##' @param max.fail integer; the maximum number of filtering failures allowed (see below).
##' If the number of filtering failures exceeds this number, execution will terminate with an error.
##' By default, \code{max.fail} is set to infinity, so no error can be triggered.
##' @param pred.mean logical; if \code{TRUE}, the prediction means are calculated for the state variables and parameters.
##' @param pred.var logical; if \code{TRUE}, the prediction variances are calculated for the state variables and parameters.
##' @param filter.mean logical; if \code{TRUE}, the filtering means are calculated for the state variables and parameters.
##' @param filter.traj logical; if \code{TRUE}, a filtered trajectory is returned for the state variables and parameters.
##' @param save.states,save.params logical.
##' If \code{save.states=TRUE}, the state-vector for each particle at each time is saved in the \code{saved.states} slot of the returned \sQuote{pfilterd_pomp} object.
##' If \code{save.params=TRUE}, the parameter-vector for each particle at each time is saved in the \code{saved.params} slot of the returned \sQuote{pfilterd_pomp} object.
##'
##' @return
##' An object of class \sQuote{pfilterd_pomp}, which extends class \sQuote{pomp}.
##'
##' @section Methods:
##' \describe{
##' \item{logLik}{ the estimated log likelihood  }
##' \item{cond.logLik}{ the estimated conditional log likelihood }
##' \item{eff.sample.size}{
##' the (time-dependent) estimated effective sample size }
##' \item{pred.mean, pred.var}{ the mean and variance of the approximate prediction distribution }
##' \item{filter.mean}{ the mean of the filtering distribution }
##' \item{filter.traj}{ retrieve one sample from the smoothing distribution}
##' \item{as.data.frame}{ coerce to a data frame }
##' \item{plot}{diagnostic plots}
##' }
##'
##' @references
##' M. S. Arulampalam, S. Maskell, N. Gordon, & T. Clapp.
##' A Tutorial on Particle Filters for Online Nonlinear, Non-Gaussian Bayesian Tracking.
##' IEEE Trans. Sig. Proc. 50:174--188, 2002.
##'
##' @examples
##'
##' pompExample(gompertz)
##' pf <- pfilter(gompertz,Np=1000)	## use 1000 particles
##' plot(pf)
##' logLik(pf)
##' cond.logLik(pf)			## conditional log-likelihoods
##' eff.sample.size(pf)             ## effective sample size
##' logLik(pfilter(pf))      	## run it again with 1000 particles
##' ## run it again with 2000 particles
##' pf <- pfilter(pf,Np=2000,filter.mean=TRUE)
##' fm <- filter.mean(pf)    	## extract the filtering means
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
    saved.params="list",
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
    saved.params=list(),
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
    stop("in ",sQuote("pfilter"),": ",sQuote("data")," is a required argument",call.=FALSE)
  }
)

setMethod(
  "pfilter",
  signature=signature(data="ANY"),
  definition=function (data, ...) {
    stop(sQuote("pfilter")," is not defined when ",sQuote("data")," is of class ",sQuote(class(data)),call.=FALSE)
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
    rinit,
    rprocess,
    dmeasure,
    params,
    Np,
    tol = 1e-17,
    max.fail = Inf,
    pred.mean = FALSE,
    pred.var = FALSE,
    filter.mean = FALSE,
    filter.traj = FALSE,
    save.states = FALSE,
    save.params = FALSE,
    ...,
    verbose = getOption("verbose", FALSE)) {

    object <- pomp(data,rinit=rinit,rprocess=rprocess,dmeasure=dmeasure,...)

    pfilter(
      object,
      params=params,
      Np=Np,
      tol=tol,
      max.fail=max.fail,
      pred.mean=pred.mean,
      pred.var=pred.var,
      filter.mean=filter.mean,
      filter.traj=filter.traj,
      save.states=save.states,
      save.params=save.params,
      verbose=verbose
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
    params,
    Np,
    tol = 1e-17,
    max.fail = Inf,
    pred.mean = FALSE,
    pred.var = FALSE,
    filter.mean = FALSE,
    filter.traj = FALSE,
    save.states = FALSE,
    save.params = FALSE,
    ...,
    verbose = getOption("verbose", FALSE)) {

    pfilter.internal(
      data,
      params=params,
      Np=Np,
      tol=tol,
      max.fail=max.fail,
      pred.mean=pred.mean,
      pred.var=pred.var,
      filter.mean=filter.mean,
      filter.traj=filter.traj,
      save.states=save.states,
      save.params=save.params,
      verbose=verbose,
      ...
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
  function (data, params, Np, tol, ...,
    verbose = getOption("verbose", FALSE)) {

    if (missing(Np)) Np <- data@Np
    if (missing(tol)) tol <- data@tol

    pfilter(as(data,"pomp"),params=params,Np=Np,tol=tol,...,
      verbose=verbose)

  }
)

pfilter.internal <- function (object, params, Np, tol, max.fail,
  pred.mean = FALSE, pred.var = FALSE, filter.mean = FALSE,
  filter.traj = FALSE, cooling, cooling.m, verbose = FALSE,
  save.states = FALSE, save.params = FALSE,
  .getnativesymbolinfo = TRUE, ...) {

  ep <- paste0("in ",sQuote("pfilter"),": ")

  object <- pomp(object,...)

  gnsi <- as.logical(.getnativesymbolinfo)
  pred.mean <- as.logical(pred.mean)
  pred.var <- as.logical(pred.var)
  filter.mean <- as.logical(filter.mean)
  filter.traj <- as.logical(filter.traj)
  verbose <- as.logical(verbose)
  tol <- as.numeric(tol)
  save.states <- as.logical(save.states)
  save.params <- as.logical(save.params)

  if (missing(params)) params <- coef(object)
  if (is.list(params)) params <- unlist(params)
  if (is.null(params)) params <- numeric(0)

  do.par.resample <- TRUE
  times <- time(object,t0=TRUE)
  ntimes <- length(times)-1

  if (missing(Np) || is.null(Np)) {
    if (is.matrix(params)) {
      Np <- ncol(params)
    } else {
      stop(ep,sQuote("Np")," must be specified.",call.=FALSE)
    }
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
  if (is.matrix(params)) {
    if (!all(Np == ncol(params)))
      stop(ep,"when ",sQuote("params")," is provided as a matrix, ",
        "you must not also specify ",sQuote("Np"),".",call.=FALSE)
  }

  Np <- as.integer(Np)

  if (length(tol) != 1 || !is.finite(tol) || tol < 0)
    stop(ep,sQuote("tol")," should be a small positive number.",call.=FALSE)

  if (NCOL(params)==1) {        # there is only one parameter vector
    do.par.resample <- FALSE
    coef(object) <- params     # set params slot to the parameters
    params <- as.matrix(params)
  }

  paramnames <- rownames(params)
  if (is.null(paramnames))
    stop(ep,sQuote("params")," must have rownames.",call.=FALSE)

  pompLoad(object,verbose=verbose)
  on.exit(pompUnload(object,verbose=verbose))

  init.x <- rinit(object,params=params,nsim=Np[1L],
    .getnativesymbolinfo=gnsi)
  statenames <- rownames(init.x)
  nvars <- nrow(init.x)
  x <- init.x

  ## set up storage for saving samples from filtering distributions
  if (save.states || filter.traj) {
    xparticles <- setNames(vector(mode="list",length=ntimes),time(object))
  }
  if (save.params) {
    pparticles <- setNames(vector(mode="list",length=ntimes),time(object))
  } else {
    pparticles <- list()
  }
  if (filter.traj) {
    pedigree <- vector(mode="list",length=ntimes+1)
  }

  loglik <- rep(NA,ntimes)
  eff.sample.size <- numeric(ntimes)
  nfail <- 0

  ## set up storage for prediction means, variances, etc.
  if (pred.mean) {
    pred.m <- matrix(
      data=0,
      nrow=nvars,
      ncol=ntimes,
      dimnames=list(
        variable=statenames,
        time=time(object))
    )
  } else {
    pred.m <- array(data=numeric(0),dim=c(0,0))
  }

  if (pred.var) {
    pred.v <- matrix(
      data=0,
      nrow=nvars,
      ncol=ntimes,
      dimnames=list(
        variable=statenames,
        time=time(object))
    )
  } else {
    pred.v <- array(data=numeric(0),dim=c(0,0))
  }

  if (filter.mean) {
    filt.m <- matrix(
      data=0,
      nrow=nvars,
      ncol=ntimes,
      dimnames=list(
        variable=statenames,
        time=time(object))
    )
  } else {
    filt.m <- array(data=numeric(0),dim=c(0,0))
  }

  if (filter.traj) {
    filt.t <- array(
      data=0,
      dim=c(nvars,1,ntimes+1),
      dimnames=list(
        variable=statenames,
        rep=1,
        time=times)
    )
  } else {
    filt.t <- array(data=numeric(0),dim=c(0,0,0))
  }

  for (nt in seq_len(ntimes)) { ## main loop

    ## advance the state variables according to the process model
    X <- tryCatch(
      rprocess(
        object,
        xstart=x,
        times=times[c(nt,nt+1)],
        params=params,
        offset=1,
        .getnativesymbolinfo=gnsi
      ),
      error = function (e) {
        stop(ep,"process simulation error: ",
          conditionMessage(e),call.=FALSE)
      }
    )

    if (pred.var) { ## check for nonfinite state variables and parameters
      problem.indices <- unique(which(!is.finite(X),arr.ind=TRUE)[,1L])
      if (length(problem.indices)>0) {  # state variables
        stop(
          ep,"non-finite state variable(s): ",
          paste(rownames(X)[problem.indices],collapse=', '),
          call.=FALSE
        )
      }
    }

    ## determine the weights
    weights <- tryCatch(
      dmeasure(
        object,
        y=object@data[,nt,drop=FALSE],
        x=X,
        times=times[nt+1],
        params=params,
        log=FALSE,
        .getnativesymbolinfo=gnsi
      ),
      error = function (e) {
        stop(ep,conditionMessage(e),call.=FALSE)
      }
    )
    if (!all(is.finite(weights))) {
      first <- which(!is.finite(weights))[1L]
      datvals <- object@data[,nt]
      weight <- weights[first]
      states <- X[,first,1L]
      params <- if (do.par.resample) params[,first] else params[,1L]
      msg <- nonfinite_dmeasure_error(time=times[nt+1],lik=weight,datvals,
        states,params)
      stop(ep,msg,call.=FALSE)
    }

    gnsi <- FALSE

    ## compute prediction mean, prediction variance, filtering mean,
    ## effective sample size, log-likelihood
    ## also do resampling if filtering has not failed
    xx <- tryCatch(
      .Call(
        pfilter_computations,
        x=X,
        params=params,
        Np=Np[nt+1],
        predmean=pred.mean,
        predvar=pred.var,
        filtmean=filter.mean,
        trackancestry=filter.traj,
        doparRS=do.par.resample,
        weights=weights,
        tol=tol
      ),
      error = function (e) {
        stop(ep,conditionMessage(e),call.=FALSE) # nocov
      }
    )
    all.fail <- xx$fail
    loglik[nt] <- xx$loglik
    eff.sample.size[nt] <- xx$ess

    x <- xx$states
    params <- xx$params

    if (pred.mean)
      pred.m[,nt] <- xx$pm
    if (pred.var)
      pred.v[,nt] <- xx$pv
    if (filter.mean)
      filt.m[,nt] <- xx$fm
    if (filter.traj)
      pedigree[[nt]] <- xx$ancestry

    if (all.fail) { ## all particles are lost
      nfail <- nfail+1
      if (verbose)
        message("filtering failure at time t = ",times[nt+1])
      if (nfail>max.fail)
        stop(ep,"too many filtering failures",call.=FALSE)
    }

    if (save.states || filter.traj) {
      xparticles[[nt]] <- x
      dimnames(xparticles[[nt]]) <- setNames(dimnames(xparticles[[nt]]),c("variable","rep"))
    }

    if (save.params) {
      pparticles[[nt]] <- params
      dimnames(pparticles[[nt]]) <- setNames(dimnames(pparticles[[nt]]),c("variable","rep"))
    }

    if (verbose && (nt%%5==0))
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
    warning(
      ep,nfail,
      ngettext(
        nfail,
        msg1=" filtering failure occurred.",
        msg2=" filtering failures occurred."
      ),
      call.=FALSE
    )

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
    saved.params=pparticles,
    Np=as.integer(Np),
    tol=tol,
    nfail=as.integer(nfail),
    loglik=sum(loglik)
  )
}

nonfinite_dmeasure_error <- function (time, lik, datvals, states, params) {
  showvals <- c(time=time,lik=lik,datvals,states,params)
  m1 <- formatC(names(showvals),preserve.width="common")
  m2 <- formatC(showvals,digits=6,width=12,format="g",
    preserve.width="common")
  paste0(
    sQuote("dmeasure")," returns non-finite value.\n",
    "likelihood, data, states, and parameters are:\n",
    paste0(m1,": ",m2,collapse="\n")
  )
}
