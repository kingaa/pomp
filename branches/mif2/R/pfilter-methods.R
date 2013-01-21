setMethod("$",signature(x="pfilterd.pomp"),function (x,name) slot(x,name))
setMethod("logLik",signature(object="pfilterd.pomp"),function(object,...)object@loglik)
setMethod("eff.sample.size",signature(object="pfilterd.pomp"),function(object,...)object@eff.sample.size)
setMethod("cond.logLik",signature(object="pfilterd.pomp"),function(object,...)object@cond.loglik)

## extract the prediction means
setMethod(
          "pred.mean",
          "pfilterd.pomp",
          function (object, pars, ...) {
            if (missing(pars)) pars <- rownames(object@pred.mean)
            object@pred.mean[pars,]
          }
          )

## extract the prediction variances
setMethod(
          "pred.var",
          "pfilterd.pomp",
          function (object, pars, ...) {
            if (missing(pars)) pars <- rownames(object@pred.var)
            object@pred.var[pars,]
          }
          )


## extract the filtering means
setMethod(
          "filter.mean",
          "pfilterd.pomp",
          function (object, pars, ...) {
            if (missing(pars)) pars <- rownames(object@filter.mean)
            object@filter.mean[pars,]
          }
          )
