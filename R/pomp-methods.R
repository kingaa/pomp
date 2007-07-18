## this file contains short methods definitions

## 'coerce' method: allows for coercion of a "pomp" object to a data-frame
setAs(
      from='pomp',
      to='data.frame',
      def = function (from) {
        x <- as.data.frame(cbind(from@times,t(from@data)))
        names(x) <- c('time',rownames(from@data))
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

