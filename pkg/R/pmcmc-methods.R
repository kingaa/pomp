## this file contains short definitions of methods for the 'pmcmc' class

## extract the estimated log likelihood
setMethod('logLik','pmcmc',function(object,...)object@loglik)

## extract the filtering means
setMethod(
          'filter.mean',
          'pmcmc',
          function (object, pars, ...) {
            if (missing(pars)) pars <- rownames(object@filter.mean)
            object@filter.mean[pars,]
          }
          )

## extract the convergence record
setMethod(
          'conv.rec',
          'pmcmc',
          function (object, pars, ...) {
            if (missing(pars)) pars <- colnames(object@conv.rec)
            object@conv.rec[,pars]
          }
          )

## plot pmcmc object
setMethod(
          "plot",
          "pmcmc",
          function (x, y = NULL, ...) {
            compare.pmcmc(x)
          }
          )

