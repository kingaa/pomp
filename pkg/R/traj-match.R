setClass(
         "traj.matched.pomp",
         contains="pomp",
         representation=representation(
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

traj.match.objfun <- function (object, params, est) {
  
  if (missing(est)) est <- character(0)
  if (!is.character(est)) stop(sQuote("est")," must be a vector of parameter names")
  if (missing(params)) params <- coef(object)
  if ((!is.numeric(params))||(is.null(names(params))))
    stop(sQuote("params")," must be a named numeric vector")
  par.est.idx <- match(est,names(params))
  if (any(is.na(par.est.idx)))
    stop("parameter(s): ",sQuote(est[is.na(par.est.idx)])," not found in ",sQuote("params"))
  params <- as.matrix(params)

  obj.fn <- function (par) {
    params[par.est.idx,] <- par
    X <- trajectory(object,params=params)
    d <- dmeasure(
                  object,
                  y=object@data,
                  x=X,
                  times=time(object),
                  params=params,
                  log=TRUE
                  )
    -sum(d)
  }

  obj.fn
}

traj.match.internal <- function (object, start, est, method, gr, eval.only, ...) {
  
  if (eval.only) {
    est <- character(0)
    guess <- numeric(0)
  } else {
    if (!is.character(est)) stop(sQuote("est")," must be a vector of parameter names")
    if (length(start)<1)
      stop(sQuote("start")," must be supplied if ",sQuote("object")," contains no parameters")
    if (!all(est%in%names(start)))
      stop(sQuote("traj.match")," error: parameters named in ",sQuote("est"),
           " must exist in ",sQuote("start"),call.=FALSE)
    guess <- start[est]
  }

  obj <- as(object,"pomp")
  coef(obj) <- start

  obj.fn <- traj.match.objfun(obj,est=est)

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
    coef(obj,est) <- unname(opt$par)
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
                    gr = NULL, eval.only = FALSE, ...) {
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
                                ...
                                )
          }
          )

setMethod(
          "traj.match",
          signature=signature(object="traj.matched.pomp"),
          function (object, start, est,
                    method = c("Nelder-Mead","subplex","SANN","BFGS","sannbox"),
                    gr = NULL, eval.only = FALSE, ...) {
            if (missing(start)) start <- coef(object)
            if (missing(est)) est <- object@est
            method <- match.arg(method)
            traj.match.internal(
                                object=object,
                                start=start,
                                est=est,
                                method=method,
                                gr=gr,
                                eval.only=eval.only,
                                ...
                                )
          }
          )
