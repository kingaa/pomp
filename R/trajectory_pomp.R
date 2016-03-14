trajectory.internal <- function (object, params, times, t0, as.data.frame = FALSE, .getnativesymbolinfo = TRUE, ...) {

  pompLoad(object)
  
  if (missing(times))
    times <- time(object,t0=FALSE)
  else
    times <- as.numeric(times)

  as.data.frame <- as.logical(as.data.frame)

  if (length(times)==0)
    stop(sQuote("times")," is empty, there is no work to do",call.=FALSE)
  
  if (any(diff(times)<=0))
    stop(sQuote("times")," must be an increasing sequence of times",call.=FALSE)

  if (missing(t0))
    t0 <- timezero(object)
  else
    t0 <- as.numeric(t0)
  
  if (t0>times[1L])
    stop("the zero-time ",sQuote("t0")," must occur no later than the first observation",call.=FALSE)
  ntimes <- length(times)
  
  if (missing(params)) {
    params <- coef(object)
    if (length(params)==0) {
      stop("trajectory error: ",sQuote("params")," must be supplied",call.=FALSE)
    }
  }
  params <- as.matrix(params)
  nrep <- ncol(params)
  paramnames <- rownames(params)
  if (is.null(paramnames))
    stop("trajectory error: ",sQuote("params")," must have rownames",call.=FALSE)

  x0 <- init.state(object,params=params,t0=t0)
  nvar <- nrow(x0)
  statenames <- rownames(x0)
  dim(x0) <- c(nvar,nrep,1)
  dimnames(x0) <- list(statenames,NULL,NULL)
  
  type <- object@skeleton.type          # map or vectorfield?
  
  if (is.na(type))
    stop("trajectory error: 'skeleton.type' unspecified",call.=FALSE)

  if (type=="map") {

    x <- .Call(iterate_map,object,times,t0,x0,params,.getnativesymbolinfo)
    .getnativesymbolinfo <- FALSE
    
  } else if (type=="vectorfield") {

    znames <- object@zeronames
    if (length(znames)>0) x0[znames,,] <- 0

    .Call(pomp_desolve_setup,object,x0,params,.getnativesymbolinfo)
    .getnativesymbolinfo <- FALSE

    X <- try(
             ode(
                 y=x0,
                 times=c(t0,times),
                 method="lsoda",
                 func="pomp_vf_eval",
                 dllname="pomp",
                 initfunc=NULL,
                 parms=NULL,
                 ...
                 ),
             silent=FALSE
             )

    .Call(pomp_desolve_takedown)

    if (inherits(X,'try-error'))
      stop("trajectory error: error in ODE integrator",call.=FALSE)
    if (attr(X,'istate')[1L]!=2)
      warning("abnormal exit from ODE integrator, istate = ",attr(X,'istate'),call.=FALSE)

    x <- array(data=t(X[-1L,-1L]),dim=c(nvar,nrep,ntimes),dimnames=list(statenames,NULL,NULL))
    
    for (z in znames)
      for (r in seq_len(ncol(x)))
        x[z,r,-1] <- diff(x[z,r,])
    
  } else {
    
    stop("deterministic skeleton not specified")

  }

  dimnames(x) <- setNames(dimnames(x),c("variable","rep","time"))

  if (as.data.frame) {
    x <- lapply(
                seq_len(ncol(x)),
                function (k) {
                  nm <- rownames(x)
                  y <- x[,k,,drop=FALSE]
                  dim(y) <- dim(y)[c(1L,3L)]
                  y <- as.data.frame(t(y))
                  names(y) <- nm
                  y$time <- times
                  y$traj <- as.integer(k)
                  y
                }
                )
    x <- do.call(rbind,x)
    x$traj <- factor(x$traj)
  }

  pompUnload(object)

  x
}

setMethod("trajectory",signature=signature(object="pomp"),
          definition=function (object, params, times, t0, as.data.frame = FALSE, ...)
          trajectory.internal(object=object,params=params,times=times,t0=t0,as.data.frame=as.data.frame,...)
          )
