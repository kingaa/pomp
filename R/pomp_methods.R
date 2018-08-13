## this file contains some basic methods definitions

## 'coerce' method: allows for coercion of a "pomp" object to a data-frame
setAs(
  from="pomp",
  to="data.frame",
  def = function (from) {
    x <- as.data.frame(cbind(from@times,t(from@data)))
    names(x) <- c(from@timename,rownames(from@data))
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
    dimnames(y) <- setNames(list(vars,time(object)),
      c("variable",object@timename))
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
      varnames <- rownames(object@states)
      if (missing(vars)) vars <- varnames
      else if (!all(vars%in%varnames))
        stop("in ",sQuote("states"),": some elements of ",
          sQuote("vars")," correspond to no state variable",call.=FALSE)
      x <- object@states[vars,,drop=FALSE]
      dimnames(x) <- setNames(list(vars,time(object)),
        c("variable",object@timename))
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
      stop(ep,sQuote("value")," must be a numeric vector.",call.=FALSE)
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
      stop(ep,"the times specified must be an increasing sequence.",call.=FALSE)
    if (object@t0>object@times[1])
      stop(ep,"the zero-time ",sQuote("t0")," must occur no later than the first observation.",call.=FALSE)
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
    if (!(is.numeric(start) && is.finite(start) && length(start)==1 &&
        is.numeric(end) && is.finite(end) && length(end)==1))
      stop("in ",sQuote("window"),": ",sQuote("start")," and ",sQuote("end"),
        " must be finite times.",call.=FALSE)
    if (!isTRUE(start <= end))
      stop("in ",sQuote("window"),": ",sQuote("start")," must not be later ",
        "than ",sQuote("end"),".",call.=FALSE)
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
    if (!(is.numeric(value) && length(value) == 1L && is.finite(value)))
      stop(ep,"the zero-time ",sQuote("t0"),
        " must be a single finite number.",call.=FALSE)
    if (value > object@times[1L])
      stop(ep,"the zero-time ",sQuote("t0"),
        " must occur no later than the first observation.",call.=FALSE)
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
        params <- partrans(object,params=object@params,dir="toEst")
      else
        params <- object@params
      if (missing(pars))
        pars <- names(params)
      else {
        excl <- setdiff(pars,names(params))
        if (length(excl)>0) {
          stop("in ",sQuote("coef"),": name(s) ",
            paste(sQuote(excl),collapse=","),
            " correspond to no parameter(s).",
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
          value <- partrans(object,params=value,dir="fromEst")
        pars <- names(value)
        if (is.null(pars)) {
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
        if (transform) val <- partrans(object,params=val,dir="fromEst")
        object@params <- val
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
          params <- partrans(object,params=params,dir="fromEst")
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
    on.exit(sink(file=NULL))
    cat("==================\npomp object ",sQuote(nm),":\n\n",sep="")
    cat("- data:\n")
    cat("  -",length(object@times),"records of",
      nrow(object@data),
      ngettext(nrow(object@data),"observable,","observables,"),
      "recorded from t =",
      min(object@times),"to",max(object@times),"\n")
    cat("  - summary of data:\n")
    print(summary(as.data.frame(t(obs(object)))))
    cat("\n- zero time, t0 = ",object@t0,"\n",sep="")
    if (length(object@tcovar)>0) {
      cat("- covariates:")
      cat("  -",nrow(object@covar),"records of",
        ncol(object@covar),"covariates,",
        "recorded from t =",min(object@tcovar),
        "to",max(object@tcovar),"\n")
      cat("  - summary of covariates:\n")
      print(summary(as.data.frame(object@covar)))
    }
    cat("- initial state simulator, rinit:\n")
    if (object@default.init) {
      cat("\t\t(default initializer)\n")
    } else {
      show(object@rinit)
    }
    cat("- process-model simulator, rprocess:\n")
    show(object@rprocess)
    cat("- process model density, dprocess:\n")
    show(object@dprocess)
    cat("- measurement model simulator, rmeasure:\n")
    show(object@rmeasure)
    cat("- measurement model density, dmeasure:\n")
    show(object@dmeasure)
    cat("- prior simulator, rprior:\n")
    show(object@rprior)
    cat("- prior density, dprior:\n")
    show(object@dprior)
    cat("- deterministic skeleton:\n")
    show(object@skeleton)
    if (object@partrans@has) {
      cat("- parameter transformations:\n")
      show(object@partrans)
    }
    if (length(coef(object))>0) {
      cat("- parameter vector:\n")
      print(coef(object))
    } else {
      cat ("- parameter vector unspecified\n");
    }
    if (length(object@userdata)>0) {
      cat("- extra user-defined variables: ",
        paste(sapply(names(object@userdata),sQuote),collapse=", "),
        "\n")
    }
    ## now display C snippets
    if (length(object@solibs) > 0) {
      for (i in seq_along(object@solibs)) {
        cat("- C snippet file ",i,":\n\n")
        cat(object@solibs[[i]]$src)
      }
    }
    file.show(f,delete.file=TRUE)
    invisible(NULL)
  }
)
