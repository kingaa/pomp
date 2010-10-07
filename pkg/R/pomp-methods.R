## this file contains some basic methods definitions

## functions to extract or call the components of a "pomp" object
setGeneric("data.array",function(object,...)standardGeneric("data.array"))

setGeneric("obs",function(object,...)standardGeneric("obsns"))

setGeneric("time<-",function(object,...,value)standardGeneric("time<-"))  

setGeneric("coef<-",function(object,...,value)standardGeneric("coef<-"))

setGeneric("states",function(object,...)standardGeneric("states"))

setGeneric("timezero",function(object,...)standardGeneric("timezero"))

setGeneric("timezero<-",function(object,...,value)standardGeneric("timezero<-"))

## 'coerce' method: allows for coercion of a "pomp" object to a data-frame
setAs(
      from="pomp",
      to="data.frame",
      def = function (from) {
        x <- as.data.frame(cbind(from@times,t(from@data)))
        names(x) <- c("time",rownames(from@data))
        if (length(from@states)>0) {
          nm <- names(x)
          x <- cbind(x,t(from@states))
          names(x) <- c(nm,rownames(from@states))
        }
        x
      }
      )

as.data.frame.pomp <- function (x, row.names, optional, ...) as(x,"data.frame")

## a simple method to extract the data array
setMethod(
          "data.array",
          "pomp",
          function (object, vars, ...) {
            varnames <- rownames(object@data)
            if (missing(vars))
              vars <- varnames
            else if (!all(vars%in%varnames))
              stop("some elements of ",sQuote("vars")," correspond to no observed variable")
            object@data[vars,,drop=FALSE]
          }
          )

## a simple method to extract the data array
setMethod(
          "obs",
          "pomp",
          function (object, vars, ...) {
            varnames <- rownames(object@data)
            if (missing(vars))
              vars <- varnames
            else if (!all(vars%in%varnames))
              stop("some elements of ",sQuote("vars")," correspond to no observed variable")
            object@data[vars,,drop=FALSE]
          }
          )

## a simple method to extract the array of states
setMethod(
          "states",
          "pomp",
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
          function (object, t0 = FALSE, ..., value) {
            if (!is.numeric(value))
              stop(sQuote("value")," must be a numeric vector",call.=TRUE)
            storage.mode(value) <- "double"
            tt0 <- object@t0
            tt <- object@times
            dd <- object@data
            ss <- object@states
            if (t0) {
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
                                     dim=c(nrow(ss),length(object@times)),
                                     dimnames=list(rownames(ss),NULL)
                                     )
              for (kt in seq_along(object@times)) {
                wr <- which(object@times[kt]==tt)
                if (length(wr)>0)
                  object@states[,kt] <- ss[,wr[1]]
              }
            }
            object
          }
          )

setMethod(
          "window",
          "pomp",
          function (x, start, end, ...) {
            tm <- time(x,t0=FALSE)
            if (missing(start))
              start <- tm[1]
            if (missing(end))
              end <- tm[length(tm)]
            tm <- tm[(tm>=start)&(tm<=end)]
            time(x,t0=FALSE) <- tm
            x
          }
          )

setMethod("timezero","pomp",function(object,...)object@t0)

setMethod(
          "timezero<-",
          "pomp",
          function(object,...,value) {
            if (value>object@times[1])
              if (!is.numeric(value) || length(value) > 1)
                stop("the zero-time ",sQuote("t0")," must be a single number",call.=TRUE)
            if (value > object@times[1])
              stop("the zero-time ",sQuote("t0")," must occur no later than the first observation",call.=TRUE)
            storage.mode(value) <- "double"
            object@t0 <- value
            object
          }
          )

## extract the coefficients
setMethod(
          "coef",
          "pomp",
          function (object, pars, ...) {
            if (missing(pars)) {
              pars <- names(object@params)
            } else {
              excl <- !(pars%in%names(object@params))
              if (any(excl)) {
                stop(
                     "in ",sQuote("coef"),": name(s) ",
                     paste(sapply(pars[excl],sQuote),collapse=","),
                     " correspond to no parameter(s)"
                     )
              }
            }
            object@params[pars]
          }
          )

## modify the coefficients
setMethod(
          "coef<-",
          "pomp",
          function (object, pars, ..., value) {
            if (missing(pars)) {          ## replace the whole params slot with 'value'
              pars <- names(value)
              if (is.null(pars))
                stop("in ",sQuote("coef<-"),": ",sQuote("value")," must be a named vector")
              object@params <- numeric(length(pars))
              names(object@params) <- pars
              object@params[] <- as.numeric(value)
            } else { ## replace or append only the parameters named in 'pars'
              if (!is.null(names(value))) ## we ignore the names of 'value'
                warning("in ",sQuote("coef<-"),": names of ",sQuote("value")," are being discarded",call.=FALSE)
              if (length(object@params)==0) { ## no pre-existing 'params' slot
                object@params <- numeric(length(pars))
                names(object@params) <- pars
                object@params[] <- as.numeric(value)
              } else { ## pre-existing params slot
                excl <- !(pars%in%names(object@params)) ## new parameters
                if (any(excl)) {                        ## append parameters
                  warning(
                          "in ",sQuote("coef<-"),": name(s) ",
                          paste(sQuote(pars[excl]),collapse=","),
                          " are not existing parameter(s);",
                          " they are being concatenated",
                          call.=FALSE
                          )
                  x <- c(object@params,numeric(length(excl)))
                  names(x) <- c(names(object@params),pars[excl])
                  object@params <- x
                }
                object@params[pars] <- as.numeric(value)
              }
            }
            object
          }
          )

setMethod(
          "print",
          "pomp",
          function (x, ...) {
            cat("data and states:\n")
            print(as(x,"data.frame"))
            cat("\ncall:\n")
            print(x@call)
            invisible(x)
          }
          )

setMethod(
          "show",
          "pomp",
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

