setClass(
         "traj.matched.pomp",
         contains="pomp",
         slots=c(
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

tmof.internal <- function (object, params, est, transform, ...) {
  
  object <- as(object,"pomp")
  if (missing(est)) est <- character(0)
  est <- as.character(est)
  transform <- as.logical(transform)

  if (missing(params)) params <- coef(object)
  if ((!is.numeric(params))||(is.null(names(params))))
    stop(sQuote("params")," must be a named numeric vector")
  if (transform)
    params <- partrans(object,params,dir="toEstimationScale")
  par.est.idx <- match(est,names(params))
  if (any(is.na(par.est.idx)))
    stop("parameter(s): ",sQuote(est[is.na(par.est.idx)])," not found in ",sQuote("params"))

  function (par) {
    pompLoad(object)
    params[par.est.idx] <- par
    if (transform)
      tparams <- partrans(object,params,dir="fromEstimationScale")
    d <- dmeasure(
                  object,
                  y=object@data,
                  x=trajectory(
                    object,
                    params=if (transform) tparams else params,
                    ...
                    ),
                  times=time(object),
                  params=if (transform) tparams else params,
                  log=TRUE
                  )
    pompUnload(object)
    -sum(d)
  }
}

setMethod(
          "traj.match.objfun",
          signature=signature(object="pomp"),
          function (object, params, est, transform = FALSE, ...)
          tmof.internal(
                        object=object,
                        params=params,
                        est=est,
                        transform=transform,
                        ...
                        )
          )

setMethod(
          "traj.match",
          signature=signature(object="pomp"),
          function (object, start, est = character(0),
                    method = c("Nelder-Mead","subplex","SANN","BFGS",
                      "sannbox","nloptr"),
                    transform = FALSE, ...)
          {

            if (missing(start)) start <- coef(object)

            method <- match.arg(method)
            est <- as.character(est)
            transform <- as.logical(transform)
            
            m <- minim.internal(
                                objfun=traj.match.objfun(
                                  object=object,
                                  params=start,
                                  est=est,
                                  transform=transform
                                  ),
                                start=start,
                                est=est,
                                object=object,
                                method=method,
                                transform=transform,
                                ...
                                )

            ## fill params slot appropriately
            coef(object) <- m$params
            
            ## fill states slot appropriately
            x <- trajectory(object)
            object@states <- array(data=x,dim=dim(x)[c(1L,3L)])
            rownames(object@states) <- rownames(x)
  
            new(
                "traj.matched.pomp",
                object,
                transform=transform,
                est=est,
                value=-m$value,
                evals=m$evals,
                convergence=m$convergence,
                msg=m$msg
                )
          }
          )

setMethod(
          "traj.match",
          signature=signature(object="traj.matched.pomp"),
          function (object, est, transform, ...)
          {
            if (missing(est)) est <- object@est
            if (missing(transform)) transform <- object@transform

            f <- selectMethod("traj.match","pomp")

            f(
              object=object,
              est=est,
              transform=transform,
              ...
              )
          }
          )
