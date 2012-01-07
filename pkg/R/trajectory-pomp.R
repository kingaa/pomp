trajectory <- function (object, params, times, t0, ...)            
  stop("function ",sQuote("trajectory")," is undefined for objects of class ",sQuote(class(object)))

setGeneric('trajectory')                                                                            

trajectory.internal <- function (object, params, times, t0, ...) {

  if (missing(times))
    times <- time(object,t0=FALSE)
  else
    times <- as.numeric(times)

  if (length(times)==0)
    stop("if ",sQuote("times")," is empty, there is no work to do",call.=FALSE)
  
  if (any(diff(times)<0))
    stop(sQuote("times")," must be a nondecreasing sequence of times",call.=FALSE)

  if (missing(t0))
    t0 <- timezero(object)
  else
    t0 <- as.numeric(t0)
  
  if (t0>times[1])
    stop("the zero-time ",sQuote("t0")," must occur no later than the first observation",call.=FALSE)
  ntimes <- length(times)
  
  if (missing(params)) {
    params <- coef(object)
    if (length(params)==0) {
      stop("trajectory error: ",sQuote("params")," must be supplied",call.=FALSE)
    }
  }
  nrep <- NCOL(params)
  if (is.null(dim(params))) {
    params <- matrix(
                     params,
                     nrow=length(params),
                     ncol=nrep,
                     dimnames=list(
                       names(params),
                       NULL
                       )
                     )
  }
  paramnames <- rownames(params)
  if (is.null(paramnames))
    stop("trajectory error: ",sQuote("params")," must have rownames",call.=FALSE)
  params <- as.matrix(params)

  x0 <- init.state(object,params=params,t0=t0)
  nvar <- nrow(x0)
  statenames <- rownames(x0)
  dim(x0) <- c(nvar,nrep,1)
  dimnames(x0) <- list(statenames,NULL,NULL)
  
  type <- object@skeleton.type          # map or vectorfield?
  
  if ("zeronames"%in%names(object@userdata))
    znames <- object@userdata$zeronames
  else
    znames <- character(0)

  if (is.na(type))
    stop("trajectory error: no skeleton specified",call.=FALSE)

  if (type=="map") {

    x <- .Call(iterate_map,object,times,t0,x0,params,znames)

  } else if (type=="vectorfield") {

    skel <- get.pomp.fun(object@skeleton)

    ## vectorfield function (RHS of ODE) in 'deSolve' format
    vf.eval <- function (t, y, ...) {
      list(
           .Call(
                 do_skeleton,
                 object,
                 x=array(
                   data=y,
                   dim=c(nvar,nrep,1),
                   dimnames=list(statenames,NULL,NULL)
                   ),
                 t=t,
                 params=params,
                 skel=skel
                 ),
           NULL
           )
    }

    if (length(znames)>0)
      x0[znames,,] <- 0

    X <- try(
             ode(
                 y=x0,
                 times=c(t0,times),
                 func=vf.eval,
                 method="lsoda",
                 ...
                 ),
             silent=FALSE
             )
    if (inherits(X,'try-error'))
      stop("trajectory error: error in ODE integrator",call.=FALSE)
    if (attr(X,'istate')[1]!=2)
      warning("abnormal exit from ODE integrator, istate = ",attr(X,'istate'),call.=FALSE)

    x <- array(data=t(X[-1,-1]),dim=c(nvar,nrep,ntimes),dimnames=list(statenames,NULL,NULL))

    if (length(znames)>0)
      x[znames,,-1] <- apply(x[znames,,,drop=FALSE],c(1,2),diff)
    
  } else {
    
    stop("deterministic skeleton not specified")

  }

  x
}

setMethod("trajectory",signature=signature(object="pomp"),definition=trajectory.internal)
