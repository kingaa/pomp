##' Forecast mean
##'
##' Mean of the one-step-ahead forecasting distribution.
##'
##' @name forecast
##' @rdname forecast
##' @aliases forecast,missing-method forecast,ANY-method
##' @family extraction methods
##' @include kalman.R
##' @inheritParams filter.mean
##'

setGeneric(
  "forecast",
  function (object, ...)
    standardGeneric("forecast")
)

setMethod(
  "forecast",
  signature=signature(object="missing"),
  definition=function (...) {
    reqd_arg("forecast","object")
  }
)

setMethod(
  "forecast",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    undef_method("forecast",object)
  }
)

##' @rdname forecast
##' @export
setMethod(
  "forecast",
  signature=signature(object="kalmand_pomp"),
  definition=function (object, vars, ...) {
    if (!missing(vars)) 
      object@forecast[vars,,drop=FALSE]
    else
      object@forecast
  }
)

##' @rdname forecast
##' @export
setMethod(
  "forecast",
  signature=signature(object="pfilterd_pomp"),
  definition=function (object, vars, ...) {
    if (undefined(object@emeasure))
      pStop("forecast",paste(sQuote(c("emeasure")),collapse=", "),
        " is a needed basic component.")
    x <- pred.mean(object)
    if (length(x)==0)
      pStop("forecast","no prediction mean. ",
        "Rerun ",sQuote("pfilter")," with ",
        sQuote("pred.mean=TRUE"),".")
    y <- emeasure(
      object,
      x=pred.mean(object),
      times=time(object),
      params=coef(object)
    )
    if (!missing(vars))
      y <- y[vars,,,drop=FALSE]
    nm <- rownames(y)
    dn <- dim(y)[c(1L,3L)]
    dim(y) <- dn
    rownames(y) <- nm
    y
  }
)
