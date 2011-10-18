## this file contains some basic methods definitions

## functions to extract or call the components of a "pomp" object
setGeneric("data.array",function(object,...)standardGeneric("data.array"))

setGeneric("obs",function(object,...)standardGeneric("obs"))

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
        if (length(from@covar)>0) {
          nm <- names(x)
          y <- .Call(lookup_in_table,from@tcovar,from@covar,from@times)
          x <- cbind(x,t(y))
          names(x) <- c(nm,rownames(y))
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
            if (length(object@states)==0) {
              NULL
            } else {
              if (missing(vars))
                vars <- seq(length=nrow(object@states))
              object@states[vars,,drop=FALSE]
            }
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

pomp.transform <- function (object, params, dir = c("forward","inverse")) {
  dir <- match.arg(dir)
  r <- length(dim(params))
  nm <- if (r>0) rownames(params) else names(params)
  tfunc <- switch(
                  dir,
                  forward=function (x) do.call(object@par.trans,c(list(x),object@userdata)),
                  inverse=function (x) do.call(object@par.untrans,c(list(x),object@userdata))
                  )
  if (r > 1) {
    retval <- apply(params,seq.int(from=2,to=r),tfunc)
    no.names <- is.null(rownames(retval))
  } else {
    retval <- tfunc(params)
    no.names <- is.null(names(retval))
  }
  if (no.names)
    switch(
           dir,
           forward=stop(
             "invalid ",sQuote("pomp")," object: ",
             sQuote("parameter.transform")," must return a named numeric vector"
             ),
           inverse=stop(
             "invalid ",sQuote("pomp")," object: ",
             sQuote("parameter.inv.transform")," must return a named numeric vector"
             )
           )
  retval
}

## extract the coefficients
setMethod(
          "coef",
          "pomp",
          function (object, pars, transform = FALSE, ...) {
            if (transform) 
              params <- pomp.transform(object,params=object@params,dir="inverse")
            else
              params <- object@params
            if (missing(pars))
              pars <- names(params)
            else {
              excl <- setdiff(pars,names(params))
              if (length(excl)>0) {
                stop(
                     "in ",sQuote("coef"),": name(s) ",
                     paste(sQuote(excl),collapse=","),
                     " correspond to no parameter(s)"
                     )
              }
            }
            params[pars]
          }
          )

## modify the coefficients
setMethod(
          "coef<-",
          "pomp",
          function (object, pars, transform = FALSE, ..., value) {
            if (missing(pars)) {          ## replace the whole params slot with 'value'
              if (transform) 
                value <- pomp.transform(object,params=value,dir="forward")
              pars <- names(value)
              if (is.null(pars)) {
                if (transform)
                  stop(sQuote("parameter.transform(value)")," must be a named vector")
                else
                  stop(sQuote("value")," must be a named vector")
              }
              object@params <- value
            } else { ## replace or append only the parameters named in 'pars'
              if (!is.null(names(value))) ## we ignore the names of 'value'
                warning(
                        "in ",sQuote("coef<-"),
                        " names of ",sQuote("value")," are being discarded",
                        call.=FALSE
                        )
##              if (length(pars)!=length(value))
##                stop(sQuote("pars")," and ",sQuote("value")," must be of equal length")
              if (length(object@params)==0) { ## no pre-existing 'params' slot
                val <- numeric(length(pars))
                names(val) <- pars
                val[] <- value
                if (transform)
                  value <- pomp.transform(object,params=val,dir="forward")
                object@params <- value
              } else { ## pre-existing params slot
                params <- coef(object,transform=transform)
                val <- numeric(length(pars))
                names(val) <- pars
                val[] <- value
                excl <- !(pars%in%names(params)) ## new parameter names
                if (any(excl)) { ## append parameters
                  warning(
                          "in ",sQuote("coef<-"),": name(s) ",
                          paste(sQuote(pars[excl]),collapse=","),
                          " do not refer to existing parameter(s);",
                          " they are being concatenated",
                          call.=FALSE
                          )
                  params <- c(params,val[excl])
                }
                params[pars] <- val
                if (transform)
                  params <- pomp.transform(object,params=params,dir="forward")
                object@params <- params
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
            cat("parameter transform function = \n")
            show(object@par.trans)
            cat("parameter inverse transform function = \n")
            show(object@par.untrans)
            if (length(object@userdata)>0) {
              cat("userdata = \n")
              show(object@userdata)
            }
            invisible(NULL)
          }
          )

