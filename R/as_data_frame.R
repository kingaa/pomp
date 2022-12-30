##' Coerce to data frame
##'
##' All \pkg{pomp} model objects can be recast as data frames.
##' The contents of the resulting data frame depend on the nature of the object.
##'
##' @name as.data.frame
##' @docType methods
##' @keywords internal
##' @rdname as_data_frame
##' @include pomp_class.R pfilter.R bsmc2.R mif2.R pmcmc.R abc.R wpfilter.R
##' @include probe.R kalman.R melt.R
##'
NULL

setAs(
  from="pomp",
  to="data.frame",
  def = function (from) {
    x <- as.data.frame(cbind(from@times,t(from@data)))
    names(x) <- c(from@timename,rownames(from@data))
    if (length(from@states) > 0) {
      nm <- c(names(x),rownames(from@states))
      x <- cbind(x,t(from@states))
      names(x) <- nm
    }
    cnames <- get_covariate_names(from@covar)
    if (length(cnames) > 0) {
      nm <- c(names(x),cnames) # see issue #56
      y <- .Call(P_lookup_in_table,from@covar,from@times)
      x <- cbind(x,t(y))
      names(x) <- nm
    }
    x
  }
)

##' @rdname as_data_frame
##' @inheritParams base::as.data.frame
##' @details
##' When \code{object} is a simple \sQuote{pomp} object,
##' \code{as(object,"data.frame")} or \code{as.data.frame(object)} results in a
##' data frame with the times, observables, states (if known), and
##' interpolated covariates (if any).
##' @export
as.data.frame.pomp <- function (x, ...) as(x,"data.frame")

setAs(
  from="pfilterd_pomp",
  to="data.frame",
  def = function (from) {
    pm <- pred_mean(from)
    pv <- pred_var(from)
    fm <- filter_mean(from)
    out <- cbind(
      as(as(from,"pomp"),"data.frame"),
      ess=eff_sample_size(from),
      cond.logLik=cond_logLik(from)
    )
    if (length(pm)>0) {
      pm <- as.data.frame(t(pm))
      names(pm) <- paste0("pred.mean.",names(pm))
      out <- cbind(out,pm)
    }
    if (length(pv)>0) {
      pv <- as.data.frame(t(pv))
      names(pv) <- paste0("pred.var.",names(pv))
      out <- cbind(out,pv)
    }
    if (length(fm)>0) {
      fm <- as.data.frame(t(fm))
      names(fm) <- paste0("filter.mean.",names(fm))
      out <- cbind(out,fm)
    }
    out
  }
)

##' @rdname as_data_frame
##' @details When \code{object} is a \sQuote{pfilterd_pomp} object,
##' coercion to a data frame results in a data frame with the same content as for a simple \sQuote{pomp},
##' but with conditional log likelihood and effective sample size estimates included, as well as filtering means, prediction means, and prediction variances, if these have been computed.
##' @export
as.data.frame.pfilterd_pomp <- function (x, ...) as(x,"data.frame")

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

##' @rdname as_data_frame
##' @details When \code{object} is a \sQuote{probed_pomp} object,
##' coercion to a data frame results in a data frame with the values of the probes computed on the data and on simulations.
##' @export
as.data.frame.probed_pomp <- function (x, ...)
  as(x,"data.frame")

setAs(
  from="kalmand_pomp",
  to="data.frame",
  def = function (from) {
    pm <- pred_mean(from)
    fm <- filter_mean(from)
    fc <- forecast(from)
    out <- cbind(
      as(as(from,"pomp"),"data.frame"),
      cond.logLik=cond_logLik(from)
    )
    if (length(pm)>0) {
      pm <- as.data.frame(t(pm))
      names(pm) <- paste0("pred.mean.",names(pm))
      out <- cbind(out,pm)
    }
    if (length(fm)>0) {
      fm <- as.data.frame(t(fm))
      names(fm) <- paste0("filter.mean.",names(fm))
      out <- cbind(out,fm)
    }
    if (length(fc)>0) {
      fc <- as.data.frame(t(fc))
      names(fc) <- paste0("forecast.",names(fc))
      out <- cbind(out,fc)
    }
    out
  }
)

##' @rdname as_data_frame
##' @details When \code{object} is a \sQuote{kalmand_pomp} object,
##' coercion to a data frame results in a data frame with prediction means, filter means and forecasts, in addition to the data.
##' @export
as.data.frame.kalmand_pomp <- function (x, ...) {
  as(x,"data.frame")
}

setAs(
  from="bsmcd_pomp",
  to="data.frame",
  def = function (from) {
    prior <- as.data.frame(t(from@prior))
    post <- as.data.frame(t(from@post))
    prior$.id <- "prior"
    post$.id <- "posterior"
    rbind(prior,post)
  }
)

##' @rdname as_data_frame
##' @details When \code{object} is a \sQuote{bsmcd_pomp} object,
##' coercion to a data frame results in a data frame with samples from the prior and posterior distribution.
##' The \code{.id} variable distinguishes them.
##' @export
as.data.frame.bsmcd_pomp <- function (x, ...) {
  as(x,"data.frame")
}

##' @importFrom dplyr bind_rows
setAs(
  from="listie",
  to="data.frame",
  def = function (from) {
    bind_rows(
      lapply(from,as,"data.frame"),
      .id=".id"
    )
  }
)

##' @rdname as_data_frame
##' @export
as.data.frame.pompList <- function (x, ...) {
  as(x,"data.frame")
}

##' @rdname as_data_frame
##' @export
as.data.frame.pfilterList <- function (x, ...) {
  as(x,"data.frame")
}

##' @rdname as_data_frame
##' @export
as.data.frame.abcList <- function (x, ...) {
  as(x,"data.frame")
}

##' @rdname as_data_frame
##' @export
as.data.frame.mif2List <- function (x, ...) {
  as(x,"data.frame")
}

##' @rdname as_data_frame
##' @export
as.data.frame.pmcmcList <- function (x, ...) {
  as(x,"data.frame")
}

setAs(
  from="wpfilterd_pomp",
  to="data.frame",
  def = function (from) {
    cbind(
      as(as(from,"pomp"),"data.frame"),
      eff.sample.size=eff_sample_size(from),
      cond.logLik=cond_logLik(from)
    )
  }
)

##' @rdname as_data_frame
##' @details When \code{object} is a \sQuote{wpfilterd_pomp} object,
##' coercion to a data frame results in a data frame with the same content as for a simple \sQuote{pomp},
##' but with conditional log likelihood and effective sample size estimates included.
##' @export
as.data.frame.wpfilterd_pomp <- function (x, ...) as(x,"data.frame")

