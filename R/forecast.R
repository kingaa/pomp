##' Forecast mean
##'
##' Mean of the one-step-ahead forecasting distribution.
##'
##' @name forecast
##' @rdname forecast
##' @aliases forecast,missing-method forecast,ANY-method
##' @family extraction methods
##' @include pfilter.R kalman.R melt.R
##' @inheritParams filter_mean
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
  definition=function (object, vars, ...,
    format = c("array", "data.frame")) {
    if (missing(vars)) {
      x <- object@forecast
    } else {
      x <- object@forecast[vars,,drop=FALSE]
    }
    format <- match.arg(format)
    if (format == "data.frame") {
      x <- melt(object@forecast[vars,,drop=FALSE])
      x$time <- time(object)[as.integer(x$time)]
    }
    x
  }
)

##' @rdname forecast
##' @export
setMethod(
  "forecast",
  signature=signature(object="pfilterd_pomp"),
  definition=function (object, vars, ...,
    format = c("array", "data.frame")) {
    if (undefined(object@emeasure))
      pStop(who="forecast",paste(sQuote(c("emeasure")),collapse=", "),
        " is a needed basic component.")
    x <- pred_mean(object)
    if (length(x)==0)
      pStop(who="forecast","no prediction mean. ",
        "Rerun ",sQuote("pfilter")," with ",
        sQuote("pred.mean=TRUE"),".")
    y <- emeasure(
      object,
      x=x,
      times=time(object),
      params=coef(object)
    )
    if (!missing(vars))
      y <- y[vars,,,drop=FALSE]
    dn <- dim(y)[c(1L,3L)]
    dnm <- dimnames(y)[c(1L,3L)]
    dim(y) <- dn
    dimnames(y) <- dnm
    format <- match.arg(format)
    if (format=="data.frame") {
      y <- melt(y)
      y$time <- time(object)[as.integer(y$time)]
    }
    y
  }
)
