##' Coerce to data frame
##'
##' All \pkg{pomp} model objects can be recast as data frames.
##' The contents of the resulting data frame depend on the nature of the object.
##'
##' @name as.data.frame
##' @rdname as_data_frame
##' @include pomp_class.R pfilter.R bsmc2.R mif2.R pmcmc.R abc.R
##' @include probe.R kalman.R
NULL

##' @name coerce-pomp-data.frame
##' @aliases coerce,pomp,data.frame-method
##' @rdname as_data_frame
##'
##' @details
##' When \code{object} is a simple \sQuote{pomp} object,
##' \code{as(object,"data.frame")} or \code{as.data.frame(object)} results in a
##' data frame with the times, observables, states (if known), and
##' interpolated covariates (if any).
##'
setAs(
  from="pomp",
  to="data.frame",
  def = function (from) {
    x <- as.data.frame(cbind(from@times,t(from@data)))
    names(x) <- c(from@timename,rownames(from@data))
    if (length(from@states)>0) {
      nm <- names(x)
      x <- cbind(x,t(from@states))
      names(x) <- c(nm,rownames(from@states))
    }
    cnames <- get_covariate_names(from@covar)
    if (length(cnames) > 0) {
      nm <- c(names(x),cnames)  # perhaps not strictly necessary (but see issue #56)
      y <- .Call(lookup_in_table,from@covar,from@times)
      x <- cbind(x,t(y))
      names(x) <- nm
    }
    x
  }
)

##' @method as.data.frame pomp
##' @rdname as_data_frame
##'
##' @param x the object to be coerced
##' @param \dots ignored
##'
as.data.frame.pomp <- function (x, ...) as(x,"data.frame")

##' @name as,pfilterd_pomp-method
##' @aliases coerce,pfilterd_pomp,data.frame-method
##' @rdname as_data_frame
##'
##' @details When \code{object} is a simple \sQuote{pfilterd_pomp} object,
##' coercion to a data frame results in a data frame with the same content as for a simple \sQuote{pomp},
##' but with conditional log likelihood and effective sample size estimates included.
##'
setAs(
  from="pfilterd_pomp",
  to="data.frame",
  def = function (from) {
    pm <- pred.mean(from)
    pv <- pred.var(from)
    fm <- filter.mean(from)
    out <- cbind(
      as(as(from,"pomp"),"data.frame"),
      ess=eff.sample.size(from),
      cond.loglik=cond.logLik(from)
    )
    if (length(pm)>0)
      out <- cbind(out,pred.mean=t(pm))
    if (length(pv)>0)
      out <- cbind(out,pred.var=t(pv))
    if (length(fm)>0)
      out <- cbind(out,filter.mean=t(fm))
    out
  }
)

##' @method as.data.frame pfilterd_pomp
##' @rdname as_data_frame
as.data.frame.pfilterd_pomp <- function (x, ...) as(x,"data.frame")

##' @name as,probed_pomp-method
##' @aliases coerce,probed_pomp,data.frame-method
##' @rdname as_data_frame
##'
##' @details When \code{object} is a \sQuote{probed_pomp} object,
##' coercion to a data frame results in a data frame with the values of the probes computed on the data and on simulations.
##'
setAs(
  from="probed_pomp",
  to="data.frame",
  def = function (from) {
    x <- as.data.frame(rbind(from@datvals,from@simvals))
    row.names(x) <- seq.int(from=0,to=nrow(x)-1)
    x$.id <- factor(c("data",rep("sim",nrow(x)-1)))
    x
  }
)

##' @method as.data.frame probed_pomp
##' @rdname as_data_frame
as.data.frame.probed_pomp <- function (x, ...)
  as(x,"data.frame")

##' @name as,kalmand_pomp-method
##' @aliases coerce,kalmand_pomp,data.frame-method
##' @rdname as_data_frame
setAs(
  from="kalmand_pomp",
  to="data.frame",
  def = function (from) {
    pm <- pred.mean(from)
    fm <- filter.mean(from)
    fc <- forecast(from)
    out <- cbind(
      as(as(from,"pomp"),"data.frame"),
      cond.loglik=cond.logLik(from)
    )
    if (length(pm)>0)
      out <- cbind(out,pred.mean=t(pm))
    if (length(fm)>0)
      out <- cbind(out,filter.mean=t(fm))
    if (length(fc)>0)
      out <- cbind(out,forecast=t(fc))
    out
  }
)

##' @method as.data.frame kalmand_pomp
##' @rdname as_data_frame
as.data.frame.kalmand_pomp <- function (x, ...) {
  as(x,"data.frame")
}
