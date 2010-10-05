trajectory <- function (object, params, times, t0, ...)
  stop("function ",sQuote("trajectory")," is undefined for objects of class ",sQuote(class(object)))
setGeneric('trajectory')

trajectory.internal <- function (object, params, times, t0, ...) {

  warn.condition <- missing(t0)
  if (warn.condition) 
    warning(
            "The default behavior of ",sQuote("trajectory")," has changed.\n",
            "See the documentation (",dQuote("pomp?trajectory"),") for details\n",
            "and set the ",sQuote("times")," and ",sQuote("t0"),
            " arguments appropriately to compensate.\n",
            call.=FALSE
            )

  if (missing(times)) {
    times <- time(object,t0=FALSE)
  } else {
    times <- as.numeric(times)
  }

  if (length(times)==0)
    stop("if ",sQuote("times")," is empty, there is no work to do",call.=FALSE)
  
  if (any(diff(times)<0))
    stop(sQuote("times")," must be a nondecreasing sequence of times",call.=FALSE)

  if (missing(t0)) {
    t0 <- timezero(object)
  } else {
    t0 <- as.numeric(t0)
  }
  
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
    stop("pfilter error: ",sQuote("params")," must have rownames",call.=FALSE)
  params <- as.matrix(params)

  tm <- t0
  x0 <- init.state(object,params=params,t0=tm)
  nm <- rownames(x0)
  dim(x0) <- c(nrow(x0),nrep,1)
  dimnames(x0) <- list(nm,NULL,NULL)
  
  x <- array(
             dim=c(nrow(x0),nrep,length(times)),
             dimnames=list(rownames(x0),NULL,NULL)
             )
  switch(
         object@skeleton.type,
         map={                # iterate the map
           for (k in seq_along(times)) {
             if (tm < times[k]) {
               while (tm < times[k]) {
                 x0[,,] <- skeleton(object,x=x0,t=tm,params=params)
                 tm <- tm+1
               }
             }
             x[,,k] <- x0
           }
         },
         vectorfield={        # integrate the vectorfield
           for (j in seq_len(nrep)) {
             X <- try(
                      deSolve::lsoda(
                                     y=x0[,j,1],
                                     times=c(t0,times),
                                     func=function(t,y,parms){
                                       list(
                                            skeleton(
                                                     object,
                                                     x=array(
                                                       data=y,
                                                       dim=c(length(y),1,1),
                                                       dimnames=list(names(y),NULL,NULL)
                                                       ),
                                                     t=t,
                                                     params=as.matrix(parms)
                                                     ),
                                            NULL
                                            )
                                     },
                                     parms=params[,j],
                                     ...
                                     ),
                      silent=FALSE
                      )
             if (inherits(X,'try-error'))
               stop("trajectory error: error in ",sQuote("lsoda"),call.=FALSE)
             if (attr(X,'istate')[[1]]!=2)
               warning("abnormal exit from ",sQuote("lsoda"),", istate = ",attr(X,'istate'),call.=FALSE)
             x[,j,] <- t(X[-1,-1])
           }
         },
         unspecified=stop("deterministic skeleton not specified")
         )
  x
}

setMethod("trajectory",signature=signature(object="pomp"),definition=trajectory.internal)
