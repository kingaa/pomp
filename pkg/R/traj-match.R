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

traj.match.internal <- function (object, start, est, method, gr, eval.only, ...) {
  
  if (eval.only) {
    par.est <- integer(0)
  } else {
    if (!is.character(est)) stop(sQuote("est")," must be a vector of parameter names")
    if (!all(est%in%names(start)))
      stop(sQuote("traj.match")," error: parameters named in ",sQuote("est"),
           " must exist in ",sQuote("start"),call.=FALSE)
    par.est <- which(names(start)%in%est)
    guess <- start[par.est]
  }

  t0 <- timezero(object)
  obj <- as(object,"pomp")

  obj.fn <- function (x) {
    p <- start
    p[par.est] <- x
    d <- dmeasure(
                  obj,
                  y=obs(obj),
                  x=trajectory(obj,params=p,t0=t0),
                  times=time(obj),
                  params=as.matrix(p),
                  log=TRUE
                  )
    -sum(d)
  }

  if (eval.only) {

    coef(obj,names(start)) <- unname(start)
    val <- obj.fn(start)
    conv <- NA
    evals <- c(1,0)
    msg <- "no optimization performed"
    
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
    msg <- if (is.null(opt$message)) character(0) else opt$message
    conv <- opt$convergence
    evals <- opt$counts
    val <- opt$value

  }

  obj@states <- trajectory(obj,t0=t0)[,1,]
  
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

setGeneric("traj.match",function(object,...)standardGeneric("traj.match"))

setMethod(
          "traj.match",
          signature=signature(object="pomp"),
          function (object, start, est,
                    method = c("Nelder-Mead","SANN","subplex"), 
                    gr = NULL, eval.only = FALSE, ...) {
            if (missing(start)) start <- coef(object)
            if (missing(est)) {
              est <- character(0)
              eval.only <- TRUE
            }
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
                    method = c("Nelder-Mead","SANN","subplex"), 
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
