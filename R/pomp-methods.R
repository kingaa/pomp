## this file contains short methods definitions

## 'coerce' method: allows for coercion of a "pomp" object to a data-frame
setAs(
      from='pomp',
      to='data.frame',
      def = function (from) {
        x <- as.data.frame(cbind(from@times,t(from@data)))
        names(x) <- c('time',rownames(from@data))
        if (length(from@states)>0) {
          nm <- names(x)
          x <- cbind(x,t(from@states[,-1,drop=FALSE]))
          names(x) <- c(nm,rownames(from@states))
        }
        x
      }
      )

## a simple method to extract the data array
setMethod(
          'data.array',
          'pomp',
          function (object, vars, ...) {
            if (missing(vars))
              vars <- seq(length=nrow(object@data))
            object@data[vars,,drop=FALSE]
          }
          )

## a simple method to extract the vector of times
setMethod('time','pomp',function(x,...)x@times)

## extract the coefficients
setMethod(
          'coef',
          'pomp',
          function (object, pars, ...) {
            if (missing(pars)) pars <- names(object@params)
            object@params[pars]
          }
          )

## modify the coefficients
setMethod(
          'coef<-',
          'pomp',
          function (object, pars, ..., value) {
            if (missing(pars)) pars <- names(value)
            if (is.null(pars)) {
              if (length(value)!=length(object@params)) {
                stop(
                     "in 'coef<-': ",
                     "number of items to replace is not equal to replacement length",
                     call.=F
                     )
              } else {
                pars <- names(value) <- names(object@params)
              }
            }
            if (is.null(names(value))) {
              if (length(value)!=length(pars)) {
                stop(
                     "in 'coef<-': ",
                     "number of items to replace is not equal to number of names given in 'pars'",
                     call.=F
                     )
              } else {
                names(value) <- pars
              }
            }
            incl <- pars%in%names(value)
            if (!all(incl)) {
              stop(
                   "in 'coef<-': ",
                   "names '",paste(pars[!incl],collapse=','),"' are not present in 'value'",
                   call.=F
                   )
            } else {
              incl <- pars%in%names(object@params)
            }
            if (!all(incl)) {
              warning(
                      "in 'coef<-': ",
                      "values '",
                      paste(pars[!incl],collapse=','),
                      "' correspond to no element in 'coef(object)' and hence are not replaced",
                      call.=F
                      )
            }
            object@params[pars[incl]] <- value[pars[incl]]
            object
          }
          )
