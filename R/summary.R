setGeneric(
    "summary",
    function (object, ...)
        standardGeneric("summary")
)

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

setMethod(
  "summary",
  "spectd_pomp",
  function (object, ...) {
    list(
      coef=coef(object),
      nsim=nrow(object@simspec),
      pvals=object@pvals
    )
  }
)

setMethod(
  "summary",
  signature=signature(object="traj_matched_pomp"),
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

setMethod(
  "summary",
  "probe_matched_pomp",
  function (object, ...) {
    c(
      summary(as(object,"probed_pomp")),
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
