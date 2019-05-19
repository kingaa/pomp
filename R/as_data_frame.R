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
      y <- .Call(P_lookup_in_table,from@covar,from@times)
      x <- cbind(x,t(y))
      names(x) <- nm
    }
    x
  }
)

##' @method as.data.frame pomp
##' @rdname as_data_frame
##'
##' @inheritParams base::as.data.frame
##' @export
##'
as.data.frame.pomp <- function (x, ...) as(x,"data.frame")

##' @name as,pfilterd_pomp-method
##' @aliases coerce,pfilterd_pomp,data.frame-method
##' @rdname as_data_frame
##'
##' @details When \code{object} is a \sQuote{pfilterd_pomp} object,
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

##' @method as.data.frame pfilterd_pomp
##' @rdname as_data_frame
##' @export
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
##' @export
as.data.frame.probed_pomp <- function (x, ...)
  as(x,"data.frame")

##' @name as,kalmand_pomp-method
##' @aliases coerce,kalmand_pomp,data.frame-method
##' @rdname as_data_frame
##'
##' @details When \code{object} is a \sQuote{kalmand_pomp} object,
##' coercion to a data frame results in a data frame with prediction means, filter means and forecasts, in addition to the data.
##'
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

##' @method as.data.frame kalmand_pomp
##' @rdname as_data_frame
##' @export
as.data.frame.kalmand_pomp <- function (x, ...) {
  as(x,"data.frame")
}

##' @name as,bsmcd_pomp-method
##' @aliases coerce,bsmcd_pomp,data.frame-method
##' @rdname as_data_frame
##'
##' @details When \code{object} is a \sQuote{bsmcd_pomp} object,
##' coercion to a data frame results in a data frame with samples from the prior and posterior distribution.
##' The \code{.id} variable distinguishes them.
##'
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

##' @method as.data.frame bsmcd_pomp
##' @rdname as_data_frame
##' @export
as.data.frame.bsmcd_pomp <- function (x, ...) {
  as(x,"data.frame")
}

##' @name as,listie-method
##' @aliases coerce,listie,data.frame-method
##' @rdname as_data_frame
##' @importFrom plyr rbind.fill
setAs(
  from="listie",
  to="data.frame",
  def = function (from) {
    x <- lapply(from,as,"data.frame")
    if (is.null(names(from))) {
      id <- seq_along(from)
    } else {
      id <- names(from)
    }
    for (k in seq_along(x)) {
      x[[k]]$.id <- id[k]
    }
    plyr::rbind.fill(x)
  }
)

##' @method as.data.frame pompList
##' @rdname as_data_frame
##' @export
as.data.frame.pompList <- function (x, ...) {
  as(x,"data.frame")
}

##' @method as.data.frame pfilterList
##' @rdname as_data_frame
##' @export
as.data.frame.pfilterList <- function (x, ...) {
  as(x,"data.frame")
}

##' @method as.data.frame abcList
##' @rdname as_data_frame
##' @export
as.data.frame.abcList <- function (x, ...) {
  as(x,"data.frame")
}

##' @method as.data.frame mif2List
##' @rdname as_data_frame
##' @export
as.data.frame.mif2List <- function (x, ...) {
  as(x,"data.frame")
}

##' @method as.data.frame pmcmcList
##' @rdname as_data_frame
##' @export
as.data.frame.pmcmcList <- function (x, ...) {
  as(x,"data.frame")
}
