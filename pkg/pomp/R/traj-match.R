setClass(
         "traj.matched.pomp",
         contains="pomp",
         representation=representation(
           start="numeric",
           transform="logical",
           est="character",
           evals="integer",
           convergence="integer",
           msg="character",
           value="numeric"
           )
         )

setMethod("$",signature=signature(x="traj.matched.pomp"),function(x, name)slot(x,name))

setMethod("logLik",signature=signature(object="traj.matched.pomp"),function(object, ...)object@value)

setMethod(
          "summary",
          signature=signature(object="traj.matched.pomp"),
          function (object, ...) {
            c(
              list(
                   params=coef(object),
                   loglik=object@value,
                   eval=object@evals,
                   convergence=object@convergence
                   ),
              if(length(object@msg)>0) list(msg=object@msg) else NULL
              )
          }
          )

traj.match.objfun <- function (object, params, est, transform = FALSE, ...) {
  
  transform <- as.logical(transform)
  if (missing(est)) est <- character(0)
  if (!is.character(est)) stop(sQuote("est")," must be a vector of parameter names")
  if (missing(params)) params <- coef(object)
  if ((!is.numeric(params))||(is.null(names(params))))
    stop(sQuote("params")," must be a named numeric vector")
  if (transform)
    params <- partrans(object,params,dir="inverse")
  par.est.idx <- match(est,names(params))
  if (any(is.na(par.est.idx)))
    stop("parameter(s): ",sQuote(est[is.na(par.est.idx)])," not found in ",sQuote("params"))

  obj.fn <- function (par) {
    params[par.est.idx] <- par
    if (transform) {
      tparams <- partrans(object,params,dir="forward")
      d <- dmeasure(
                    object,
                    y=object@data,
                    x=trajectory(object,params=tparams,...),
                    times=time(object),
                    params=tparams,
                    log=TRUE
                    )
    } else {
      d <- dmeasure(
                    object,
                    y=object@data,
                    x=trajectory(object,params=params,...),
                    times=time(object),
                    params=params,
                    log=TRUE
                    )
    }
    -sum(d)
  }

  obj.fn
}

traj.match.internal <- function (object, start, est, method, gr, eval.only, transform, ...) {
  
  transform <- as.logical(transform)

  if (eval.only) {
    est <- character(0)
    guess <- numeric(0)
    transform <- FALSE
  } else {
    if (!is.character(est)) stop(sQuote("est")," must be a vector of parameter names")
    if (length(start)<1)
      stop(sQuote("start")," must be supplied if ",sQuote("object")," contains no parameters")
    if (transform) {
      tstart <- partrans(object,start,dir="inverse")
      if (is.null(names(tstart))||(!all(est%in%names(tstart))))
        stop(sQuote("est")," must refer to parameters named in ",sQuote("partrans(object,start,dir=\"inverse\")"))
      guess <- tstart[est]
    } else {
      if (is.null(names(start))||(!all(est%in%names(start))))
        stop(sQuote("est")," must refer to parameters named in ",sQuote("start"))
      guess <- start[est]
    }
  }

  obj <- as(object,"pomp")
  coef(obj) <- start

  obj.fn <- traj.match.objfun(obj,est=est,transform=transform)

  if (eval.only) {

    val <- obj.fn(guess)
    conv <- NA
    evals <- c(1,0)
    msg <- "no optimization performed"
    
  } else {

    if (method=="subplex") {
      opt <- subplex::subplex(par=guess,fn=obj.fn,control=list(...))
    } else if (method=="sannbox") {
      opt <- sannbox(par=guess,fn=obj.fn,control=list(...))
    } else {
      opt <- optim(par=guess,fn=obj.fn,gr=gr,method=method,control=list(...))
    }

    if (!is.null(names(opt$par)) && !all(est==names(opt$par)))
      stop("mismatch between parameter names returned by optimizer and ",sQuote("est"))
    coef(obj,est,transform=transform) <- unname(opt$par)
    msg <- if (is.null(opt$message)) character(0) else opt$message
    conv <- opt$convergence
    evals <- opt$counts
    val <- opt$value

  }

  ## fill 'states' slot of returned object with the trajectory
  x <- trajectory(obj)
  obj@states <- array(data=x,dim=dim(x)[c(1,3)])
  rownames(obj@states) <- rownames(x)
  
  new(
      "traj.matched.pomp",
      obj,
      start=start,
      transform=transform,
      est=as.character(est),
      evals=as.integer(evals),
      convergence=as.integer(conv),
      msg=msg,
      value=as.numeric(-val)
      )
}

traj.match <- function (object, ...)
  stop("function ",sQuote("traj.match")," is undefined for objects of class ",sQuote(class(object)))

setGeneric("traj.match")

setMethod(
          "traj.match",
          signature=signature(object="pomp"),
          function (object, start, est,
                    method = c("Nelder-Mead","subplex","SANN","BFGS","sannbox"),
                    gr = NULL, eval.only = FALSE, transform = FALSE, ...) {
            transform <- as.logical(transform)
            if (missing(start)) start <- coef(object)
            if (!eval.only && missing(est))
              stop(sQuote("est")," must be supplied if optimization is to be done")
            if (eval.only) est <- character(0)
            method <- match.arg(method)
            traj.match.internal(
                                object=object,
                                start=start,
                                est=est,
                                method=method,
                                gr=gr,
                                eval.only=eval.only,
                                transform=transform,
                                ...
                                )
          }
          )

setMethod(
          "traj.match",
          signature=signature(object="traj.matched.pomp"),
          function (object, start, est,
                    method = c("Nelder-Mead","subplex","SANN","BFGS","sannbox"),
                    gr = NULL, eval.only = FALSE, transform, ...) {
            if (missing(start)) start <- coef(object)
            if (missing(est)) est <- object@est
            if (missing(transform)) transform <- object@transform
            transform <- as.logical(transform)
            method <- match.arg(method)
            traj.match.internal(
                                object=object,
                                start=start,
                                est=est,
                                method=method,
                                gr=gr,
                                eval.only=eval.only,
                                transform=transform,
                                ...
                                )
          }
          )
