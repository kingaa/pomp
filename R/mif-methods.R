## this file contains short definitions of methods for the 'mif' class

## extract the estimated log likelihood
setMethod('logLik','mif',function(object,...)object@loglik)

## extract the coefficients
setMethod(
          'coef',
          'mif',
          function (object, pars, ...) {
            if (missing(pars)) pars <- names(object@coef)
            object@coef[pars]
          }
          )

## modify the coefficients
setMethod(
          'coef<-',
          'mif',
          function (object, pars, ..., value) {
            if (missing(pars)) pars <- names(value)
            if (is.null(pars)) {
              if (length(value)!=length(object@coef))
                stop(
                     "in 'coef<-': ",
                     "number of items to replace is not equal to replacement length",
                     call.=F
                     )
              pars <- names(value) <- names(object@coef)
            }
            if (is.null(names(value))) {
              if (length(value)!=length(pars))
                stop(
                     "in 'coef<-': ",
                     "number of items to replace is not equal to number of names given in 'pars'",
                     call.=F
                     )
              names(value) <- pars
            }
            incl <- pars%in%names(value)
            if (!all(incl))
              stop(
                   "in 'coef<-': ",
                   "names '",paste(pars[!incl],collapse=','),"' are not present in 'value'",
                   call.=F
                   )
            incl <- pars%in%names(object@coef)
            if (!all(incl))
              warning(
                      "in 'coef<-': ",
                      "values '",
                      paste(pars[!incl],collapse=','),
                      "' correspond to no element in 'coef(object)' and hence are not replaced",
                      call.=F
                      )
            object@coef[pars[incl]] <- value[pars[incl]]
            object
          }
          )

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
