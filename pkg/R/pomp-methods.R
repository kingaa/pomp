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

## a simple method to extract the array of states
setMethod(
          'states',
          'pomp',
          function (object, vars, ...) {
            if (missing(vars))
              vars <- seq(length=nrow(object@states))
            object@states[vars,,drop=FALSE]
          }
          )

## a simple method to extract the vector of times
setMethod(
          "time",
          "pomp",
          function (x, t0 = FALSE, ...) {
            if (t0) c(x@t0,x@times) else x@times
          }
          )

## extract the coefficients
setMethod(
          'coef',
          'pomp',
          function (object, pars, ...) {
            if (missing(pars)) {
              pars <- names(object@params)
            } else {
              excl <- !(pars%in%names(object@params))
              if (any(excl)) {
                stop(
                     "in ",sQuote("coef"),": name(s) ",
                     paste(sapply(pars[excl],sQuote),collapse=','),
                     " correspond to no parameter(s)"
                     )
              }
            }
            object@params[pars]
          }
          )

## modify the coefficients
setMethod(
          'coef<-',
          'pomp',
          function (object, pars, ..., value) {
            if (length(object@params)==0) {
              if (missing(pars)) {
                pars <- names(value)
                if (is.null(pars))
                  stop("in ",sQuote("coef<-"),": ",sQuote("value")," must be a named vector")
              } else {
                if (length(pars)!=length(value))
                  stop("in ",sQuote("coef<-"),": ",sQuote("pars")," and ",sQuote("value")," must be of the same length")
              }
              object@params <- as.numeric(value)
              names(object@params) <- pars
            } else {
              if (missing(pars)) {
                pars <- names(object@params)
                if (is.null(pars))
                  stop("bad ",sQuote("pomp")," object: slot ",sQuote("params")," should be a named vector")
              } else {
                excl <- !(pars%in%names(object@params))
                if (any(excl)) {
                  stop(
                       "in ",sQuote("coef<-"),": name(s) ",
                       paste(sapply(pars[excl],sQuote),collapse=','),
                       " correspond to no parameter(s)"
                       )
                }
              }
              if (length(pars)!=length(value))
                stop("in ",sQuote("coef<-"),": ",sQuote("pars")," and ",sQuote("value")," must be of the same length")
              object@params[pars] <- as.numeric(value)
            }
            object
          }
          )

setMethod(
          'print',
          'pomp',
          function (x, ...) {
            print(as(x,'data.frame'))
            invisible(x)
          }
          )

setMethod(
          'show',
          'pomp',
          function (object) {
            print(object)
            cat("zero time, t0 = ",object@t0,"\n")
            if (length(coef(object))>0) {
              cat("parameter(s):\n")
              print(coef(object))
            } else {
              cat ("parameter(s) unspecified\n");
            }
            cat("process model simulator, rprocess = \n")
            print(object@rprocess)
            cat("process model density, dprocess = \n")
            print(object@dprocess)
            cat("measurement model simulator, rmeasure = \n")
            print(object@rmeasure)
            cat("measurement model density, dmeasure = \n")
            print(object@dmeasure)
            cat("initializer = \n")
            print(object@initializer)
            cat("userdata = \n")
            show(object@userdata)
            invisible(NULL)
          }
          )

