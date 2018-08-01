## initialize the state variables of the process model

setMethod(
  "init.state",
  signature=signature("pomp"),
  definition=function (object, params, t0, nsim, ...) {
    init.state.internal(object=object,params=params,t0=t0,nsim=nsim,...)
  }
)

init.state.internal <- function (object, params, t0, nsim,
  .getnativesymbolinfo = TRUE, ...) {
  if (missing(t0)) t0 <- object@t0
  if (missing(params)) params <- coef(object)
  else storage.mode(params) <- "double"
  if (missing(nsim)) nsim <- NCOL(params)
  pompLoad(object)
  x <- .Call(do_init_state,object,params,t0,nsim,.getnativesymbolinfo)
  pompUnload(object)
  x
}
