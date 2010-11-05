setClass(
         "traj.matched.pomp",
         contains="pomp",
         representation=representation(
           evals="integer",
           convergence="integer",
           msg="character",
           value="numeric"
           )
         )

setMethod("$",signature(x="traj.matched.pomp"),function(x, name) slot(x,name))

setMethod("logLik",signature(object="traj.matched.pomp"),function(object,...)object@value)

setMethod(
          "summary",
          "traj.matched.pomp",
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

trajectory.nll <- function (par, est, object, params, t0, fail.value = NA, ...) {
  if (missing(par)) par <- numeric(0)
  if (missing(est)) est <- integer(0)
  if (missing(params)) params <- coef(object)
  if (missing(t0)) t0 <- timezero(object)
  params[est] <- par
  x <- trajectory(object,params=params,t0=t0)
  d <- dmeasure(
                object,
                y=data.array(object),
                x=x,
                times=time(object),
                params=as.matrix(params),
                log=TRUE
                )
  -sum(d)
}

trajectory.ls <- function (par, est, object, params, t0, fail.value = NA, transform, ...) {
  if (missing(transform)) transform <- identity
  if (missing(par)) par <- numeric(0)
  if (missing(est)) est <- integer(0)
  if (missing(params)) params <- coef(object)
  if (missing(t0)) t0 <- timezero(object)
  params[est] <- par
  x <- trajectory(object,params=params,t0=t0)
  d <- dmeasure(
                object,
                y=data.array(object),
                x=x,
                times=time(object),
                params=as.matrix(params),
                log=TRUE
                )
  -sum(d)
}

traj.match <- function (object, start, est, method = c("Nelder-Mead","SANN","subplex"), 
                        gr = NULL, eval.only = FALSE, ...) {
  
  if (!is(object,'pomp'))
    stop("traj.match error: ",sQuote("object")," must be a ",sQuote("pomp")," object",call.=FALSE)

  if (missing(start)) start <- coef(object)
  
  method <- match.arg(method)

  if (!is.character(est)) stop(sQuote("est")," must be a vector of parameter names")
  if (!all(est%in%names(start)))
    stop("traj.match error: parameters named in ",sQuote("est")," must exist in ",sQuote("start"),call.=FALSE)
  par.est <- which(names(start)%in%est)

  guess <- start[par.est]
  t0 <- timezero(object)
  obj <- as(object,"pomp")

  obj.fn <- function (x) {
    p <- start
    p[par.est] <- x
    x <- trajectory(obj,params=p,t0=t0)
    d <- dmeasure(
                  obj,
                  y=data.array(object),
                  x=trajectory(obj,params=p,t0=t0),
                  times=time(object),
                  params=as.matrix(p),
                  log=TRUE
                  )
    -sum(d)
  }

  if (eval.only) {

    coef(obj,names(start)) <- unname(start)
    val <- obj.fn(guess)
    conv <- NA
    evals <- c(1,0)
    msg <- paste("no optimization performed")
    
  } else {

    if (method=="subplex") {

      opt <- subplex::subplex(
                              par=guess,
                              fn=obj.fn,
                              control=list(...)
                              )

    } else {

      opt <- optim(
                   par=guess,
                   fn=obj.fn,
                   gr=gr,
                   method=method,
                   control=list(...)
                   )
    }

    coef(obj,names(opt$par)) <- unname(opt$par)
    obj@states[,] <- trajectory(obj,t0=t0)[,1,]
    msg <- if (is.null(opt$message)) character(0) else opt$message
    conv <- opt$convergence
    evals <- opt$counts
    val <- opt$value

  }

  ans <- new(
             "traj.matched.pomp",
             obj,
             evals=as.integer(evals),
             convergence=as.integer(conv),
             msg=msg,
             value=as.numeric(-val)
	     )
}
