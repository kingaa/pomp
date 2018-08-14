##' Summary methods
##'
##' Display a summary of a fitted model object.
##'
##' @rdname summary
##' @name summary
##' @include probe_match.R probe.R spect.R spect_match.R traj_match.R
##'
##' @aliases summary,probed_pomp-method summary,spectd_pomp-method
##' summary,spect_matched_pomp-method summary,traj_matched_pomp-method
##'
##' @param object a fitted model object
##' @param \dots ignored
NULL

setGeneric(
    "summary",
    function (object, ...)
        standardGeneric("summary")
)

##' @name summary-probed_pomp
##' @rdname summary
setMethod(
  "summary",
  signature=signature(object="probed_pomp"),
  definition=function (object, ...) {
    list(
      coef=coef(object),
      nsim=nrow(object@simvals),
      quantiles=object@quantiles,
      pvals=object@pvals,
      synth.loglik=object@synth.loglik
    )
  }
)

##' @name summary-spectd_pomp
##' @rdname summary
setMethod(
  "summary",
  signature=signature(object="spectd_pomp"),
  definition=function (object, ...) {
    list(
      coef=coef(object),
      nsim=nrow(object@simspec),
      pvals=object@pvals
    )
  }
)

##' @name summary-traj_matched_pomp
##' @rdname summary
setMethod(
  "summary",
  signature=signature(object="traj_matched_pomp"),
  definition=function (object, ...) {
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

##' @name summary-spect_matched_pomp
##' @rdname summary
setMethod(
  "summary",
  signature=signature(object="spect_matched_pomp"),
  definition=function (object, ...) {

    c(
      summary(as(object,"spectd_pomp")),
      list(
        est=object@est,
        value=object@value,
        eval=object@evals,
        convergence=object@convergence
      ),
      if(length(object@msg)>0) list(msg=object@msg) else NULL
    )

  }
)
