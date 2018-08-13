## initialize the state variables of the process model

setMethod(
  "rinit",
  signature=signature("pomp"),
  definition=function (object, params, t0, nsim, ...) {
    rinit.internal(object=object,params=params,t0=t0,nsim=nsim,...)
  }
)

rinit.internal <- function (object, params, t0, nsim,
  .getnativesymbolinfo = TRUE, ...) {
  if (missing(t0)) t0 <- object@t0
  if (missing(params)) params <- coef(object)
  else storage.mode(params) <- "double"
  if (missing(nsim)) nsim <- NCOL(params)
  pompLoad(object)
  on.exit(pompUnload(object))
  x <- .Call(do_rinit,object,params,t0,nsim,.getnativesymbolinfo)
  x
}
