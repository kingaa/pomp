## this file contains some basic methods definitions

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
      nm <- c(names(x),colnames(from@covar))
      y <- .Call(lookup_in_table,from@tcovar,from@covar,from@times)
      x <- cbind(x,t(y))
      names(x) <- nm
    }
    x
  }
)

as.data.frame.pomp <- function (x, row.names, optional, ...) as(x,"data.frame")

## parameter transformations
setMethod(
  "partrans",
  signature=signature(object="pomp"),
  definition=function (object, params,
    dir = c("fromEstimationScale", "toEstimationScale"), ...) {
    dir <- match.arg(dir)
    partrans.internal(object=object,params=params,dir=dir,...)
  }
)

partrans.internal <- function (object, params,
  dir = c("fromEstimationScale", "toEstimationScale"),
  .getnativesymbolinfo = TRUE, ...) {
  if (!object@has.trans) return(params)
  dir <- switch(dir,fromEstimationScale=1L,toEstimationScale=-1L)
  pompLoad(object)
  rv <- .Call(do_partrans,object,params,dir,.getnativesymbolinfo)
  pompUnload(object)
  rv
}

## a simple method to extract the data array
setMethod(
  "obs",
  signature=signature(object="pomp"),
  definition=function (object, vars, ...) {
    varnames <- rownames(object@data)
    if (missing(vars))
      vars <- varnames
    else if (!all(vars%in%varnames))
      stop("in ",sQuote("obs"),": some elements of ",
        sQuote("vars")," correspond to no observed variable",call.=FALSE)
    y <- object@data[vars,,drop=FALSE]
    dimnames(y) <- list(variable=rownames(y),time=time(object))
    y
  }
)

## a simple method to extract the array of states
setMethod(
  "states",
  signature=signature(object="pomp"),
  definition=function (object, vars, ...) {
    if (length(object@states)==0) {
      NULL
    } else {
      if (missing(vars))
        vars <- seq(length=nrow(object@states))
      x <- object@states[vars,,drop=FALSE]
      dimnames(x) <- list(variable=rownames(x),time=time(object))
      x
    }
  }
)

## a simple method to extract the vector of times
setMethod(
  "time",
  signature=signature(x="pomp"),
  definition=function (x, t0 = FALSE, ...) {
    if (t0) c(x@t0,x@times) else x@times
  }
)

## modify the time
setMethod(
  "time<-",
  signature=signature(object="pomp"),
  definition=function (object, t0 = FALSE, ..., value) {
    ep <- paste0("in ",sQuote("time<-"),": ")
    if (!is.numeric(value))
      stop(ep,sQuote("value")," must be a numeric vector",call.=FALSE)
    storage.mode(value) <- "double"
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
      stop(ep,"the times specified must be an increasing sequence",call.=FALSE)
    if (object@t0>object@times[1])
      stop(ep,"the zero-time ",sQuote("t0")," must occur no later than the first observation",call.=FALSE)
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
  signature=signature(x="pomp"),
  definition=function (x, start, end, ...) {
    tm <- time(x,t0=FALSE)
    if (missing(start))
      start <- tm[1L]
    if (missing(end))
      end <- tm[length(tm)]
    tm <- tm[(tm>=start)&(tm<=end)]
    time(x,t0=FALSE) <- tm
    x
  }
)

setMethod(
  "timezero",
  signature=signature(object="pomp"),
  definition=function(object,...)object@t0
)

setMethod(
  "timezero<-",
  signature=signature(object="pomp"),
  definition=function(object,...,value) {
    ep <- paste0("in ",sQuote("timezero<-"),": ")
    if (value[1L]>object@times[1L])
      if (!is.numeric(value) || length(value) > 1L)
        stop(ep,"the zero-time ",sQuote("t0"),
          " must be a single number",call.=FALSE)
    if (value > object@times[1L])
      stop(ep,"the zero-time ",sQuote("t0"),
        " must occur no later than the first observation",call.=FALSE)
    storage.mode(value) <- "double"
    object@t0 <- value
    object
  }
)

## extract the coefficients
setMethod(
  "coef",
  signature=signature(object="pomp"),
  definition=function (object, pars, transform = FALSE, ...) {
    if (length(object@params)>0) {
      if (transform)
        params <- partrans(object,params=object@params,dir="toEstimationScale")
      else
        params <- object@params
      if (missing(pars))
        pars <- names(params)
      else {
        excl <- setdiff(pars,names(params))
        if (length(excl)>0) {
          stop("in ",sQuote("coef"),": name(s) ",
            paste(sQuote(excl),collapse=","),
            " correspond to no parameter(s)",
            call.=FALSE)
        }
      }
      params[pars]
    } else {
      numeric(0)
    }
  }
)

## modify the coefficients
setMethod(
  "coef<-",
  signature=signature(object="pomp"),
  definition=function (object, pars, transform = FALSE, ..., value) {
    ep <- paste0("in ",sQuote("coef<-"),": ")
    if (is.null(value)) value <- numeric(0)
    if (is.list(value)) value <- unlist(value)
    if (missing(pars)) {          ## replace the whole params slot with 'value'
      if (length(value)>0) {
        if (transform)
          value <- partrans(object,params=value,dir="fromEstimationScale")
        pars <- names(value)
        if (is.null(pars)) {
          if (transform)
            stop(ep,"parameter.transform(",sQuote("value"),
              ") must be a named vector",call.=FALSE)
          else
            stop(ep,sQuote("value")," must be a named vector",call.=FALSE)
        }
      }
      object@params <- value
    } else { ## replace or append only the parameters named in 'pars'
      if (!is.null(names(value))) ## we ignore the names of 'value'
        warning(ep," names of ",sQuote("value")," are being discarded",call.=FALSE)
      if (length(object@params)==0) { ## no pre-existing 'params' slot
        val <- numeric(length(pars))
        names(val) <- pars
        val[] <- value
        if (transform)
          value <- partrans(object,params=val,dir="fromEstimationScale")
        object@params <- value
      } else { ## pre-existing params slot
        params <- coef(object,transform=transform)
        val <- numeric(length(pars))
        names(val) <- pars
        val[] <- value
        excl <- !(pars%in%names(params)) ## new parameter names
        if (any(excl)) { ## append parameters
          warning(ep,"name(s) ",
            paste(sQuote(pars[excl]),collapse=","),
            " do not refer to existing parameter(s);",
            " they are being concatenated",
            call.=FALSE)
          params <- c(params,val[excl])
        }
        params[pars] <- val
        if (transform)
          params <- partrans(object,params=params,dir="fromEstimationScale")
        object@params <- params
      }
    }
    storage.mode(object@params) <- "double"
    object
  }
)

setMethod(
  "spy",
  signature=signature(object="pomp"),
  definition=function (object) {
    nm <- deparse(substitute(object,env=parent.frame()))
    f <- tempfile()
    sink(file=f)
    cat("==================\npomp object ",sQuote(nm),":\n\n",sep="")
    cat("-",length(object@times),"records of",
      nrow(obs(object)),
      ngettext(nrow(obs(object)),"observable,","observables,"),
      "recorded from t =",
      min(object@times),"to",max(object@times),"\n")
    cat("- zero time, t0 = ",object@t0,"\n",sep="")
    cat("- summary of data:\n")
    print(summary(as.data.frame(t(obs(object)))))
    if (length(object@tcovar)>0) {
      cat("-",nrow(object@covar),"records of",
        ncol(object@covar),"covariates,",
        "recorded from t =",min(object@tcovar),
        "to",max(object@tcovar),"\n")
      cat("- summary of covariates:\n")
      print(summary(as.data.frame(object@covar)))
    }
    if (object@rprocess@type == 0L) {
      cat("- undefined process-model simulator\n")
    } else if (object@rprocess@type == 4L) {
      cat("- Gillespie-method process-model simulator, rate.fun = ")
      show(object@rprocess@rate.fn)
    } else {
      if (object@rprocess@type == 1L) {
        cat("- one-step process-model simulator, step.fun = ")
      } else if (object@rprocess@type == 2L) {
        cat("- discrete-time process-model simulator, step.fun = ")
      } else {
        cat("- Euler-method process-model simulator, step.fun = ")
      }
      show(object@rprocess@step.fn)
    }
    cat("- process model density, dprocess = ")
    show(object@dprocess)
    cat("- measurement model simulator, rmeasure = ")
    show(object@rmeasure)
    cat("- measurement model density, dmeasure = ")
    show(object@dmeasure)
    cat("- prior simulator, rprior = ")
    show(object@rprior)
    cat("- prior density, dprior = ")
    show(object@dprior)
    cat("- skeleton ",
      if (object@skeleton.type!="undef")
        paste0("(",object@skeleton.type,") ")
      else "",
      "= ",sep="")
    show(object@skeleton)
    cat("- initializer = ")
    show(object@initializer)
    cat("- parameter transformation (to estimation scale) = ")
    show(object@to.trans)
    cat("- parameter transformation (from estimation scale) = ")
    show(object@from.trans)
    if (length(coef(object))>0) {
      cat("- parameter(s):\n")
      print(coef(object))
    } else {
      cat ("- parameter(s) unspecified\n");
    }
    if (length(object@userdata)>0) {
      cat("- extra user-defined variables: ",
        paste(sapply(names(object@userdata),sQuote),collapse=", "),
        "\n")
    }
    ## now display C snippets
    if (length(object@solibs) > 0) {
      cat("- C snippet file 1:\n\n")
      for (i in seq_along(object@solibs)) {
        cat(object@solibs[[i]]$src)
      }
    }
    sink(file=NULL)
    file.show(f)
    file.remove(f)
    invisible(NULL)
  }
)

## useful for prepending the 'lib' generated by 'hitch' to the 'solibs'
## slot of a 'pomp' object
setMethod(
  "solibs<-",
  signature=signature(object="pomp"),
  definition=function (object, ..., value) {
      object@solibs <- c(list(value),object@solibs)
      object
  }
)
