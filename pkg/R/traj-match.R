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
  coef(obj,names(start)) <- unname(start)
  pmat <- as.matrix(start)

  obj.fn <- function (x, object, params, t0, ind) {
    params[ind,] <- x
    X <- trajectory(object,params=params,t0=t0)
    d <- dmeasure(
                  object,
                  y=obs(object),
                  x=X,
                  times=time(object),
                  params=params,
                  log=TRUE
                  )
    -sum(d)
  }

  if (eval.only) {

    val <- obj.fn(numeric(0),object=obj,params=pmat,t0=t0,ind=par.est)
    conv <- NA
    evals <- c(1,0)
    msg <- "no optimization performed"
    
  } else {

    if (method=="subplex") {

      opt <- subplex::subplex(
                              par=guess,
                              fn=obj.fn,
                              control=list(...),
                              object=obj,
                              params=pmat,
                              t0=t0,
                              ind=par.est
                              )

    } else if (method=="sannbox") {

      opt <- sannbox(
                     par=guess,
                     fn=obj.fn,
                     control=list(...),
                     object=obj,
                     params=pmat,
                     t0=t0,
                     ind=par.est
                     )

    } else {

      opt <- optim(
                   par=guess,
                   fn=obj.fn,
                   gr=gr,
                   method=method,
                   control=list(...),
                   object=obj,
                   params=pmat,
                   t0=t0,
                   ind=par.est
                   )
      
    }

    coef(obj,names(opt$par)) <- unname(opt$par)
    msg <- if (is.null(opt$message)) character(0) else opt$message
    conv <- opt$convergence
    evals <- opt$counts
    val <- opt$value

  }

  ## fill 'states' slot of returned object with the trajectory
  x <- trajectory(obj,t0=t0)
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
                    method = c("Nelder-Mead","sannbox","subplex"), 
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
                    method = c("Nelder-Mead","sannbox","subplex"), 
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
