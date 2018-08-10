## particle filtering codes

setClass(
  "pfilterd.pomp",
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

setMethod(
  "pfilter",
  signature=signature(object="pomp"),
  function (object, params, Np,
    tol = 1e-17,
    max.fail = Inf,
    pred.mean = FALSE,
    pred.var = FALSE,
    filter.mean = FALSE,
    filter.traj = FALSE,
    save.states = FALSE,
    save.params = FALSE,
    verbose = getOption("verbose"),
    ...) {
    if (missing(params)) params <- coef(object)
    pfilter.internal(
      object=object,
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

setMethod(
  "pfilter",
  signature=signature(object="pfilterd.pomp"),
  function (object, params, Np, tol, ...) {
    if (missing(params)) params <- coef(object)
    if (missing(Np)) Np <- object@Np
    if (missing(tol)) tol <- object@tol
    f <- selectMethod("pfilter","pomp")
    f(
      object=object,
      params=params,
      Np=Np,
      tol=tol,
      ...
    )
  }
)

pfilter.internal <- function (object, params, Np,
  tol, max.fail,
  pred.mean = FALSE,
  pred.var = FALSE,
  filter.mean = FALSE,
  filter.traj = FALSE,
  cooling, cooling.m,
  verbose = FALSE,
  save.states = FALSE,
  save.params = FALSE,
  .getnativesymbolinfo = TRUE) {

  ep <- paste0("in ",sQuote("pfilter"),": ")

  object <- as(object,"pomp")
  pompLoad(object,verbose=verbose)
  on.exit(pompUnload(object,verbose=verbose))

  gnsi.rproc <- gnsi.dmeas <- as.logical(.getnativesymbolinfo)
  pred.mean <- as.logical(pred.mean)
  pred.var <- as.logical(pred.var)
  filter.mean <- as.logical(filter.mean)
  filter.traj <- as.logical(filter.traj)
  verbose <- as.logical(verbose)
  save.states <- as.logical(save.states)
  save.params <- as.logical(save.params)

  if (is.list(params)) params <- unlist(params)
  if (length(params)==0)
    stop(ep,sQuote("params")," must be specified",call.=FALSE)

  one.par <- FALSE
  times <- time(object,t0=TRUE)
  ntimes <- length(times)-1

  if (missing(Np)) {
    if (is.matrix(params)) {
      Np <- ncol(params)
    } else {
      stop(ep,sQuote("Np")," must be specified",call.=FALSE)
    }
  }
  if (is.function(Np)) {
    Np <- tryCatch(
      vapply(seq.int(from=0,to=ntimes,by=1),Np,numeric(1)),
      error = function (e) {
        stop(ep,"if ",sQuote("Np")," is a function, ",
          "it must return a single positive integer",call.=FALSE)
      }
    )
  }
  if (length(Np)==1)
    Np <- rep(Np,times=ntimes+1)
  else if (length(Np)!=(ntimes+1))
    stop(ep,sQuote("Np")," must have length 1 or length ",ntimes+1,call.=FALSE)
  if (any(Np<=0))
    stop(ep,"number of particles, ",sQuote("Np"),", must always be positive",call.=FALSE)
  if (!is.numeric(Np))
    stop(ep,sQuote("Np")," must be a number, a vector of numbers, or a function",call.=FALSE)
  Np <- as.integer(Np)
  if (is.matrix(params)) {
    if (!all(Np==ncol(params)))
      stop(ep,"when ",sQuote("params")," is provided as a matrix, do not specify ",
        sQuote("Np"),"!",call.=FALSE)
  }

  if (NCOL(params)==1) {        # there is only one parameter vector
    one.par <- TRUE
    coef(object) <- params     # set params slot to the parameters
    params <- as.matrix(params)
  }

  paramnames <- rownames(params)
  if (is.null(paramnames))
    stop(ep,sQuote("params")," must have rownames",call.=FALSE)

  init.x <- init.state(object,params=params,nsim=Np[1L])
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
        .getnativesymbolinfo=gnsi.rproc
      ),
      error = function (e) {
        stop(ep,"process simulation error: ",
          conditionMessage(e),call.=FALSE)
      }
    )
    gnsi.rproc <- FALSE

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
        .getnativesymbolinfo=gnsi.dmeas
      ),
      error = function (e) {
        stop(ep,"error in calculation of weights: ",
          conditionMessage(e),call.=FALSE)
      }
    )
    if (!all(is.finite(weights))) {
      first <- which(!is.finite(weights))[1L]
      datvals <- object@data[,nt]
      weight <- weights[first]
      states <- X[,first,1L]
      params <- if (one.par) params[,1L] else params[,first]
      msg <- nonfinite_dmeasure_error(time=times[nt+1],lik=weight,datvals,states,params)
      stop(ep,msg,call.=FALSE)
    }
    gnsi.dmeas <- FALSE

    ## compute prediction mean, prediction variance, filtering mean,
    ## effective sample size, log-likelihood
    ## also do resampling if filtering has not failed
    xx <- tryCatch(
      .Call(
        pfilter_computations,
        x=X,
        params=params,
        Np=Np[nt+1],
        rw_sd=numeric(0),
        predmean=pred.mean,
        predvar=pred.var,
        filtmean=filter.mean,
        trackancestry=filter.traj,
        onepar=one.par,
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
    "pfilterd.pomp",
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
