## this file contains some basic methods definitions

## functions to extract or call the components of a "pomp" object
data.array <- function (object, ...)
  stop("function ",sQuote("data.array")," is undefined for objects of class ",sQuote(class(object)))
setGeneric('data.array')  

"time<-" <- function (object, ..., value)
  stop("function ",sQuote("time<-")," is undefined for objects of class ",sQuote(class(object)))
setGeneric("time<-")  

"coef<-" <- function (object, pars, ..., value)
  stop("function ",sQuote("coef<-")," is undefined for objects of class ",sQuote(class(object)))
setGeneric("coef<-")

states <- function (object, ...)
  stop("function ",sQuote("states")," is undefined for objects of class ",sQuote(class(object)))
setGeneric('states')

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

## modify the time
setMethod(
          "time<-",
          "pomp",
          function (object, include.t0 = FALSE, ..., value) {
            if (!is.numeric(value))
              stop(sQuote("value")," must be a numeric vector",call.=TRUE)
            storage.mode(value) <- "double"
            tt0 <- object@t0
            tt <- object@times
            dd <- object@data
            ss <- object@states
            if (include.t0) {
              object@t0 <- value[1]
              object@times <- value[-1]
            } else {
              object@times <- value
            }
            if (!all(diff(object@times)>0))
              stop("the times specified must be an increasing sequence",call.=TRUE)
            if (object@t0>object@times[1])
              stop("the zero-time ",sQuote("t0")," must occur no later than the first observation",call.=TRUE)
            object@data <- array(
                                 data=NA,
                                 dim=c(nrow(dd),length(object@times)),
                                 dimnames=list(rownames(dd),NULL)
                                 )
            object@data[,object@times%in%tt] <- dd[,tt%in%object@times]
            if (length(ss)>0) {
               object@states <- array(
                                     data=NA,
                                     dim=c(nrow(ss),length(object@times)+1),
                                     dimnames=list(rownames(ss),NULL)
                                     )
               if (ncol(ss)>length(tt)) tt <- c(tt0,tt)
               nt <- c(object@t0,object@times)
               for (kt in seq_len(length(nt))) {
                 wr <- which(nt[kt]==tt)
                 if (length(wr)>0)
                   object@states[,kt] <- ss[,wr[1]]
               }
             }
            object
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
            cat("data and states:\n")
            print(as(x,'data.frame'))
            cat("\ncall:\n")
            print(x@call)
            invisible(x)
          }
          )

setMethod(
          'show',
          'pomp',
          function (object) {
            print(object)
            cat("zero time, t0 = ",object@t0,"\n",sep="")
            if (length(coef(object))>0) {
              cat("parameter(s):\n")
              print(coef(object))
            } else {
              cat ("parameter(s) unspecified\n");
            }
            cat("process model simulator, rprocess = \n")
            show(object@rprocess)
            cat("process model density, dprocess = \n")
            show(object@dprocess)
            cat("measurement model simulator, rmeasure = \n")
            show(object@rmeasure)
            cat("measurement model density, dmeasure = \n")
            show(object@dmeasure)
            if (!is.na(object@skeleton.type)) {
              cat("skeleton (",object@skeleton.type,") = \n")
              show(object@skeleton)
            }
            cat("initializer = \n")
            show(object@initializer)
            if (length(object@userdata)>0) {
              cat("userdata = \n")
              show(object@userdata)
            }
            invisible(NULL)
          }
          )

