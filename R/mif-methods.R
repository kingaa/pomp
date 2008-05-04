## this file contains short definitions of methods for the 'mif' class

## extract the estimated log likelihood
setMethod('logLik','mif',function(object,...)object@loglik)

## extract the prediction means
setMethod(
          'pred.mean',
          'mif',
          function (object, pars, ...) {
            if (missing(pars)) pars <- rownames(object@pred.mean)
            object@pred.mean[pars,]
          }
          )

## extract the prediction variances
setMethod(
          'pred.var',
          'mif',
          function (object, pars, ...) {
            if (missing(pars)) pars <- rownames(object@pred.var)
            object@pred.var[pars,]
          }
          )


## extract the filtering means
setMethod(
          'filter.mean',
          'mif',
          function (object, pars, ...) {
            if (missing(pars)) pars <- rownames(object@filter.mean)
            object@filter.mean[pars,]
          }
          )

## extract the convergence record
setMethod(
          'conv.rec',
          'mif',
          function (object, pars, ...) {
            if (missing(pars)) pars <- colnames(object@conv.rec)
            object@conv.rec[,pars]
          }
          )
