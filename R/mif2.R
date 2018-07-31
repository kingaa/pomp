## IF2 algorithm functions

## define the mif2d.pomp class
setClass(
  'mif2d.pomp',
  contains='pfilterd.pomp',
  slots=c(
    Nmif = 'integer',
    rw.sd = 'matrix',
    cooling.type = 'character',
    cooling.fraction.50 = 'numeric',
    transform = 'logical',
    traces = 'matrix'
  )
)

setMethod(
  "mif2",
  signature=signature(object="pomp"),
  definition = function (object, Nmif = 1, start, Np, rw.sd, transform = FALSE,
    cooling.type = c("hyperbolic", "geometric"), cooling.fraction.50,
    tol = 1e-17, max.fail = Inf, verbose = getOption("verbose"),...) {

    mif2.internal(
      object=object,
      Nmif=Nmif,
      start=start,
      Np=Np,
      rw.sd=rw.sd,
      transform=transform,
      cooling.type=cooling.type,
      cooling.fraction.50=cooling.fraction.50,
      tol=tol,
      max.fail=max.fail,
      verbose=verbose,
      ...
    )

  }
)

setMethod(
  "mif2",
  signature=signature(object="pfilterd.pomp"),
  definition = function (object, Nmif = 1, Np, tol, ...) {

    if (missing(Np)) Np <- object@Np
    if (missing(tol)) tol <- object@tol

    mif2(as(object,"pomp"),Nmif=Nmif,Np=Np,tol=tol,...)
  }
)

setMethod(
  "mif2",
  signature=signature(object="mif2d.pomp"),
  definition = function (object, Nmif, start, Np, rw.sd, transform,
    cooling.type, cooling.fraction.50, tol, ...) {

    if (missing(Nmif)) Nmif <- object@Nmif
    if (missing(start)) start <- coef(object)
    if (missing(rw.sd)) rw.sd <- object@rw.sd
    if (missing(transform)) transform <- object@transform
    if (missing(cooling.type)) cooling.type <- object@cooling.type
    if (missing(cooling.fraction.50)) cooling.fraction.50 <- object@cooling.fraction.50

    if (missing(Np)) Np <- object@Np
    if (missing(tol)) tol <- object@tol

    mif2(as(object,"pomp"),Nmif=Nmif,start=start,Np=Np,rw.sd=rw.sd,
      transform=transform,cooling.type=cooling.type,
      cooling.fraction.50=cooling.fraction.50,tol=tol,...)
  }
)

setMethod(
  'continue',
  signature=signature(object='mif2d.pomp'),
  definition = function (object, Nmif = 1, ...) {

    ndone <- object@Nmif

    obj <- mif2(object,Nmif=Nmif,...,
      .ndone=ndone,.paramMatrix=object@paramMatrix)

    object@traces[ndone+1,c('loglik','nfail')] <- obj@traces[1L,c('loglik','nfail')]
    obj@traces <- rbind(
      object@traces,
      obj@traces[-1L,colnames(object@traces)]
    )
    names(dimnames(obj@traces)) <- c("iteration","variable")
    obj@Nmif <- as.integer(ndone+Nmif)

    obj
  }
)

mif2.internal <- function (object, Nmif, start, Np, rw.sd, transform = FALSE,
  cooling.type = c("hyperbolic", "geometric"), cooling.fraction.50,
  tol = 1e-17, max.fail = Inf, verbose = FALSE, .ndone = 0L,
  .indices = integer(0), .paramMatrix = NULL,
  .getnativesymbolinfo = TRUE, ...) {

  ep <- paste0("in ",sQuote("mif2"),": ")

  transform <- as.logical(transform)
  verbose <- as.logical(verbose)
  gnsi <- as.logical(.getnativesymbolinfo)

  if (length(Nmif) != 1 || !is.numeric(Nmif) || !is.finite(Nmif) || Nmif < 1)
    stop(ep,sQuote("Nmif")," must be a positive integer.",call.=FALSE)
  Nmif <- as.integer(Nmif)

  if (missing(start)) start <- coef(object)
  if (is.list(start)) start <- unlist(start)
  if (length(start)==0 || !is.numeric(start) || is.null(names(start)))
    stop(ep,"parameters must be specified as a named numeric vector.",call.=FALSE)

  ntimes <- length(time(object))

  if (is.null(Np)) {
    stop(ep,sQuote("Np")," must be specified.",call.=FALSE)
  } else if (is.function(Np)) {
    Np <- tryCatch(
      vapply(seq_len(ntimes),Np,numeric(1)),
      error = function (e) {
        stop(ep,"if ",sQuote("Np"),
          " is a function, it must return a single positive integer.",
          call.=FALSE)
      }
    )
  } else if (!is.numeric(Np)) {
    stop(ep,sQuote("Np"),
      " must be a number, a vector of numbers, or a function.",call.=FALSE)
  }

  if (length(Np) == 1) {
    Np <- rep(Np,times=ntimes)
  } else if (length(Np) > ntimes) {
    if (Np[1L] != Np[ntimes+1] || length(Np) > ntimes+1) {
      warning(ep,"Np[k] ignored for k > ",sQuote("length(time(object))"),".",
        call.=FALSE)
    }
    Np <- head(Np,ntimes)
  } else if (length(Np) < ntimes) {
    stop(ep,sQuote("Np")," must have length 1 or ",
      sQuote("length(time(object))"),".",call.=FALSE)
  }

  if (!all(is.finite(Np)) || any(Np <= 0))
    stop(ep,"number of particles, ",sQuote("Np"),
      ", must be a positive integer.",call.=FALSE)

  Np <- as.integer(Np)
  Np <- c(Np,Np[1L])

  if (missing(rw.sd))
    stop(ep,sQuote("rw.sd")," must be specified!",call.=FALSE)
  rw.sd <- pkern.sd(rw.sd,time=time(object),paramnames=names(start))

  if (missing(cooling.fraction.50))
    stop(ep,sQuote("cooling.fraction.50")," is a required argument.",call.=FALSE)
  if (length(cooling.fraction.50) != 1 || !is.numeric(cooling.fraction.50) ||
      !is.finite(cooling.fraction.50) || cooling.fraction.50 <= 0 ||
      cooling.fraction.50 > 1)
    stop(ep,sQuote("cooling.fraction.50")," must be in (0,1]",call.=FALSE)
  cooling.fraction.50 <- as.numeric(cooling.fraction.50)
  cooling.type <- match.arg(cooling.type)

  cooling.fn <- mif2.cooling(
    type=cooling.type,
    fraction=cooling.fraction.50,
    ntimes=length(time(object))
  )

  if (is.null(.paramMatrix)) {
    paramMatrix <- array(data=start,dim=c(length(start),Np[1L]),
      dimnames=list(variable=names(start),rep=NULL))
  } else {
    paramMatrix <- .paramMatrix
    start <- apply(paramMatrix,1L,mean)
  }

  traces <- array(dim=c(Nmif+1,length(start)+2),
    dimnames=list(iteration=seq.int(.ndone,.ndone+Nmif),
      variable=c('loglik','nfail',names(start))))
  traces[1L,] <- c(NA,NA,start)

  object <- pomp(object,...)

  pompLoad(object,verbose=verbose)
  on.exit(pompUnload(object,verbose=verbose))

  if (transform)
    paramMatrix <- partrans(object,paramMatrix,dir="toEstimationScale",
      .getnativesymbolinfo=gnsi)

  ## iterate the filtering
  for (n in seq_len(Nmif)) {

    pfp <- tryCatch(
      mif2.pfilter(
        object=object,
        params=paramMatrix,
        Np=Np,
        mifiter=.ndone+n,
        cooling.fn=cooling.fn,
        rw.sd=rw.sd,
        tol=tol,
        max.fail=max.fail,
        verbose=verbose,
        transform=transform,
        .indices=.indices,
        .getnativesymbolinfo=gnsi
      ),
      error = function (e) {
        stop(ep,conditionMessage(e),call.=FALSE)
      }
    )

    gnsi <- FALSE

    paramMatrix <- pfp@paramMatrix
    traces[n+1,-c(1,2)] <- coef(pfp)
    traces[n,c(1,2)] <- c(pfp@loglik,pfp@nfail)
    .indices <- pfp@indices

    if (verbose) cat("mif2 iteration",n,"of",Nmif,"completed\n")

  }

  if (transform)
    pfp@paramMatrix <- partrans(object,paramMatrix,dir="fromEstimationScale",
      .getnativesymbolinfo=gnsi)

  new(
    "mif2d.pomp",
    pfp,
    Nmif=Nmif,
    rw.sd=rw.sd,
    cooling.type=cooling.type,
    cooling.fraction.50=cooling.fraction.50,
    transform=transform,
    traces=traces
  )
}

mif2.cooling <- function (type, fraction, ntimes) {
  switch(
    type,
    geometric={
      factor <- fraction^(1/50)
      function (nt, m) {
        alpha <- factor^(nt/ntimes+m-1)
        list(alpha=alpha,gamma=alpha^2)
      }
    },
    hyperbolic={
      if (fraction < 1) {
        scal <- (50*ntimes*fraction-1)/(1-fraction)
        function (nt, m) {
          alpha <- (1+scal)/(scal+nt+ntimes*(m-1))
          list(alpha=alpha,gamma=alpha^2)
        }
      } else {
        function (nt, m) {
          list(alpha=1,gamma=1)
        }
      }
    }
  )
}

mif2.pfilter <- function (object, params, Np, mifiter, rw.sd, cooling.fn,
  tol = 1e-17, max.fail = Inf, transform, verbose, .indices = integer(0),
  .getnativesymbolinfo = TRUE) {

  ep <- paste0("in ",sQuote("mif2.pfilter"),": ")

  tol <- as.numeric(tol)
  gnsi <- as.logical(.getnativesymbolinfo)
  transform <- as.logical(transform)
  verbose <- as.logical(verbose)
  mifiter <- as.integer(mifiter)
  Np <- as.integer(Np)

  if (length(tol) != 1 || !is.finite(tol) || tol < 0)
    stop(ep,sQuote("tol")," should be a small positive number.",call.=FALSE)

  do_ta <- length(.indices)>0L
  if (do_ta && length(.indices)!=Np[1L])
    stop(ep,sQuote(".indices")," has improper length",call.=FALSE)

  times <- time(object,t0=TRUE)
  ntimes <- length(times)-1

  loglik <- rep(NA,ntimes)
  eff.sample.size <- numeric(ntimes)
  nfail <- 0

  for (nt in seq_len(ntimes)) {

    ## perturb parameters
    pmag <- cooling.fn(nt,mifiter)$alpha*rw.sd[,nt]
    params <- .Call(randwalk_perturbation,params,pmag)

    if (transform)
      tparams <- partrans(object,params,dir="fromEstimationScale",
        .getnativesymbolinfo=gnsi)

    if (nt == 1L) {
      ## get initial states
      x <- init.state(object,params=if (transform) tparams else params)
    }

    ## advance the state variables according to the process model
    X <- tryCatch(
      rprocess(
        object,
        xstart=x,
        times=times[c(nt,nt+1)],
        params=if (transform) tparams else params,
        offset=1,
        .getnativesymbolinfo=gnsi
      ),
      error = function (e) {
        stop(ep,"process simulation error: ",
          conditionMessage(e),call.=FALSE)
      }
    )

    ## determine the weights
    weights <- tryCatch(
      dmeasure(
        object,
        y=object@data[,nt,drop=FALSE],
        x=X,
        times=times[nt+1],
        params=if (transform) tparams else params,
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
      params <- if (transform) tparams[,first] else params[,first]
      msg <- nonfinite_dmeasure_error(time=times[nt+1],lik=weight,datvals,states,params)
      stop(ep,msg,call.=FALSE)
    }
    gnsi <- FALSE

    ## compute weighted mean at last timestep
    if (nt == ntimes) {
      if (any(weights>0)) {
        coef(object,transform=transform) <- apply(params,1L,weighted.mean,w=weights)
      } else {
        warning(ep,"filtering failure at last filter iteration; using ",
          "unweighted mean for point estimate.",call.=FALSE)
        coef(object,transform=transform) <- apply(params,1L,mean)
      }
    }

    ## compute effective sample size, log-likelihood
    ## also do resampling if filtering has not failed
    xx <- tryCatch(
      .Call(
        pfilter_computations,
        x=X,
        params=params,
        Np=Np[nt+1],
        predmean=FALSE,
        predvar=FALSE,
        filtmean=FALSE,
        trackancestry=do_ta,
        doparRS=TRUE,
        weights=weights,
        tol=tol
      ),
      error = function (e) {
        stop(ep,"particle-filter error: ",conditionMessage(e),call.=FALSE) # nocov
      }
    )
    all.fail <- xx$fail
    loglik[nt] <- xx$loglik
    eff.sample.size[nt] <- xx$ess
    if (do_ta) {
      .indices <- .indices[xx$ancestry]
    }

    x <- xx$states
    params <- xx$params

    if (all.fail) { ## all particles are lost
      nfail <- nfail+1
      if (verbose)
        message("filtering failure at time t = ",times[nt+1])
      if (nfail>max.fail)
        stop(ep,"too many filtering failures",call.=FALSE)
    }

    if (verbose && (nt%%5==0))
      cat("mif2 pfilter timestep",nt,"of",ntimes,"finished\n")

  }

  if (nfail>0) {
    warning(
      ep,nfail,
      ngettext(
        nfail,
        msg1=" filtering failure occurred.",
        msg2=" filtering failures occurred."
      ),
      call.=FALSE
    )
  }

  new(
    "pfilterd.pomp",
    as(object,"pomp"),
    paramMatrix=params,
    eff.sample.size=eff.sample.size,
    cond.loglik=loglik,
    indices=.indices,
    Np=Np,
    tol=tol,
    nfail=as.integer(nfail),
    loglik=sum(loglik)
  )
}

rw.sd <- safecall

pkern.sd <- function (rw.sd, time, paramnames) {
  ep <- paste0("in ",sQuote("mif2"),": ")
  if (is.matrix(rw.sd)) return(rw.sd)
  if (is(rw.sd,"safecall")) {
    enclos <- rw.sd@envir
    rw.sd <- as.list(rw.sd@call)[-1L]
  } else {
    stop(ep,sQuote("rw.sd")," should be specified using the ",sQuote("rw.sd"),
      " function. See ",sQuote("?mif2"),".",call.=FALSE)
  }
  if (is.null(names(rw.sd)) | any(names(rw.sd)==""))
    stop(ep,"in ",sQuote("rw.sd"),": parameters must be referenced by name.",call.=FALSE)
  if (!all(names(rw.sd) %in% paramnames)) {
    unrec <- names(rw.sd)[!names(rw.sd) %in% paramnames]
    stop(ep,"the following parameter(s), ",
      "given random walks in ",sQuote("rw.sd"),", are not present in ",
      sQuote("start"),": ",paste(sapply(unrec,sQuote),collapse=","),
      call.=FALSE)
  }
  ivp <- function (sd, lag = 1L) {
    sd*(seq_along(time)==lag)
  }
  sds <- lapply(rw.sd,eval,envir=list(time=time,ivp=ivp),enclos=enclos)
  for (n in names(sds)) {
    len <- length(sds[[n]])
    if (len==1) {
      sds[[n]] <- rep(sds[[n]],length(time))
    } else if (len!=length(time)) {
      stop(ep,sQuote("rw.sd")," spec for parameter ",sQuote(n),
        " does not evaluate to a vector of the correct length (",
        sQuote("length(time(object))"),"=",length(time),").",call.=FALSE)
    }
  }
  do.call(rbind,sds)
}
