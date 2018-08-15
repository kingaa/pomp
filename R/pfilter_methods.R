## methods for the 'pfilterd_pomp' class

setMethod(
  "logLik",
  signature=signature(object="pfilterd_pomp"),
  definition=function(object,...)object@loglik
)

setMethod(
  "eff.sample.size",
  signature=signature(object="pfilterd_pomp"),
  definition=function(object,...)object@eff.sample.size
)

setMethod(
  "cond.logLik",
  signature=signature(object="pfilterd_pomp"),
  definition=function(object,...)object@cond.loglik
)

## 'coerce' method: allows for coercion of a 'pomp' object to a data-frame
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

as.data.frame.pfilterd_pomp <- function (x, row.names, optional, ...) as(x,"data.frame")

## extract the prediction means
setMethod(
  "pred.mean",
  signature=signature(object="pfilterd_pomp"),
  definition=function (object, vars, ...) {
    if (missing(vars)) vars <- rownames(object@pred.mean)
    object@pred.mean[vars,,drop=FALSE]
  }
)

## extract the prediction variances
setMethod(
  "pred.var",
  signature=signature(object="pfilterd_pomp"),
  definition=function (object, vars, ...) {
    if (missing(vars)) vars <- rownames(object@pred.var)
    object@pred.var[vars,,drop=FALSE]
  }
)

## extract the filtering means
setMethod(
  "filter.mean",
  signature=signature(object="pfilterd_pomp"),
  definition=function (object, vars, ...) {
    if (missing(vars)) vars <- rownames(object@filter.mean)
    object@filter.mean[vars,,drop=FALSE]
  }
)

setClass(
  "pfilterList",
  contains="list",
  validity=function (object) {
    if (length(object) > 0) {
      if (!all(vapply(object,is,logical(1),"pfilterd_pomp"))) {
        retval <- paste0(
          "error in ",sQuote("c"),
          ": dissimilar objects cannot be combined"
        )
        return(retval)
      }
    }
    TRUE
  }
)

setClassUnion("Pfilter",c("pfilterd_pomp","pfilterList"))

setMethod(
  "concat",
  signature=signature(...="Pfilter"),
  definition=function (...) {
    y <- lapply(
      list(...),
      function (z) {
        if (is(z,"list"))
          setNames(as(z,"list"),names(z))
        else
          z
      }
    )
    new("pfilterList",unlist(y))
  }
)

c.Pfilter <- concat

## extract the filtered trajectory
setMethod(
  "filter.traj",
  signature=signature(object="pfilterd_pomp"),
  definition=function (object, vars, ...) {
    if (missing(vars)) vars <- rownames(object@filter.traj)
    object@filter.traj[vars,,,drop=FALSE]
  }
)

setMethod(
  "filter.traj",
  signature=signature(object="pfilterList"),
  definition=function (object, ...) {
    fts <- lapply(object,filter.traj,...)
    d <- dim(fts[[1]])
    nm <- dimnames(fts[[1]])
    x <- do.call(c,fts)
    dim(x) <- c(d,length(fts))
    dimnames(x) <- c(nm,list(chain=names(fts)))
    x
  }
)

