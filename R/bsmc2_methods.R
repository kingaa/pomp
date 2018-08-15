setMethod(
  "logEvidence",
  signature=signature(object="bsmcd_pomp"),
  definition=function(object,...)object@log.evidence
)

setMethod(
  "cond.logEvidence",
  signature=signature(object="bsmcd_pomp"),
  definition=function(object,...)object@cond.log.evidence
)

setMethod(
  "prior_samples",
  signature=signature(object="bsmcd_pomp"),
  definition=function (object, pars, ...) {
    if (missing(pars)) pars <- object@est
    object@prior[pars,,drop=FALSE]
  })

setMethod(
  "posterior_samples",
  signature=signature(object="bsmcd_pomp"),
  definition=function (object, pars, ...) {
    if (missing(pars)) pars <- object@est
    object@post[pars,,drop=FALSE]
  })
