##' Constructor of the basic pomp object
##'
##' This function constructs a \sQuote{pomp} object, encoding a
##' partially-observed Markov process model together with a uni- or
##' multi-variate time series.  As such, it is central to all the package's
##' functionality.  One implements the model by specifying some or all of its
##' \emph{basic components}.  These include:
##' \describe{
##' \item{rprocess,}{the
##' simulator of the unobserved Markov state process;}
##' \item{dprocess,}{the
##' evaluator of the probability density function for transitions of the
##' unobserved Markov state process;}
##' \item{rmeasure,}{the simulator of the
##' observed process, conditional on the unobserved state;}
##' \item{dmeasure,}{the evaluator of the measurement model probability density
##' function;}
##' \item{rinit,}{which samples from the distribution of the state
##' process at the zero-time;}
##' \item{rprior,}{which samples from a prior
##' probability distribution on the parameters;}
##' \item{dprior}{which evaluates
##' the prior probability density function;}
##' \item{skeleton}{which computes the
##' deterministic skeleton of the unobserved state process.} }
##' The basic structure and its rationale are described in the \emph{Journal of
##' Statistical Software} paper, an updated version of which is to be found on
##' the \href{https://kingaa.github.io/pomp/}{package website}.
##'
##' @name pomp
##' @rdname pomp
##' @include pomp_class.R pomp_fun.R csnippet.R safecall.R builder.R
##'
##' @param data,times required; the time series data and times at which
##' observations are made.  \code{data} should be given as a data-frame and
##' \code{times} must indicate the column of observation times by name or
##' index.  \code{times} must be numeric and strictly increasing.  Internally,
##' \code{data} will be internally coerced to an array with storage-mode
##' \code{double}.
##'
##' In addition, a \sQuote{pomp} object can be supplied in the \code{data}
##' argument.  In this case, the call to \code{pomp} will add element to, or
##' replace elements of, the supplied \sQuote{pomp} object.
##' @param t0 The zero-time, at which the stochastic dynamical system is to be
##' initialized.  This must be no later than the time of the first observation,
##' i.e., \code{t0 <= times[1]}.  This argument is required whenever
##' \code{data} is a data-frame.
##' @param rinit optional; draws from the distribution of initial values of the
##' unobserved Markov state process.  Specifically, given a vector of
##' parameters, \code{params} and an initial time, \code{t0}, the rinit
##' determines the state vector at time \code{t0}.  See below under \dQuote{The
##' State-Process Initializer} for details.
##' @param rprocess,dprocess optional; specification of the simulator and
##' probability density evaluation function of the unobserved state process.
##' See below under \dQuote{The Unobserved Markov State-Process Model} for
##' details.
##'
##' \strong{Note:} it is not typically necessary (or even feasible) to define
##' \code{dprocess}.  In fact, no current \pkg{pomp} inference algorithm makes
##' use of \code{dprocess}.  This functionality is provided only to support
##' future algorithm development.
##' @param rmeasure,dmeasure optional; specifications of the measurement model.
##' See below under \dQuote{The Measurement Model} for details.
##' @param skeleton optional; the deterministic skeleton of the unobserved
##' state process.  See below under \dQuote{The Deterministic Skeleton} for
##' details.
##' @param rprior,dprior optional; specification of the prior distribution on
##' parameters.  See below under \dQuote{Specifying a Prior} for details.
##' @param partrans optional parameter transformations.  Many algorithms for
##' parameter estimation search an unconstrained space of parameters.  When
##' working with such an algorithm and a model for which the parameters are
##' constrained, it can be useful to transform parameters.  One should supply
##' the \code{partrans} argument via a call to
##' \code{parameter_trans(toEst,fromEst,\dots{},log,logit,barycentric)}.  See
##' below under \dQuote{Parameter Transformations} for more details.
##' @param params optional; named numeric vector of parameters.  This will be
##' coerced internally to storage mode \code{double}.
##' @param covar optional covariate table, constructed using
##' \code{covariate_table}.
##'
##' If a covariate table is supplied, then the value of each of the covariates
##' is interpolated as needed.  The resulting interpolated values are made
##' available to the appropriate basic components.  Note that \code{covar} will
##' be coerced internally to storage mode \code{double}.  See below under
##' \dQuote{Covariates} for more details.
##' @param obsnames,statenames,paramnames,covarnames optional character vectors
##' specifying the names of observables, state variables, parameters, and
##' covariates, respectively.  These are used only in the event that one or
##' more of the basic components are defined using C snippets or native
##' routines.  It is usually unnecessary to specify \code{obsnames} or
##' \code{covarnames}, as these will by default be read from \code{data} and
##' \code{covars}, respectively.
##' @param zeronames optional character vector specifying the names of
##' accumulator variables (see below under \dQuote{Accumulator Variables}).
##' @param PACKAGE optional string giving the name of the dynamically loaded
##' library in which any native routines are to be found.  This is only useful
##' if one or more of the model components has been specified using a
##' precompiled dynamically loaded library; it is not used for any component
##' specified using C snippets.
##' @param globals optional character; C code that will be included in the
##' source for (and therefore hard-coded into) the shared-object library
##' created when the call to \code{pomp} uses C snippets.  If no C snippets are
##' used, \code{globals} has no effect.
##' @param cdir,cfile,shlib.args optional character variables.  \code{cdir}
##' specifies the name of the directory within which C snippet code will be
##' compiled.  By default, this is in a temporary directory specific to the
##' running instance of .  \code{cfile} gives the name of the file (in
##' directory \code{cdir}) into which C snippet codes will be written.  By
##' default, a random filename is used.  The \code{shlib.args} can be used to
##' pass command-line arguments to the \code{R CMD SHLIB} call that will
##' compile the C snippets.
##' @param \dots Any additional arguments given to \code{pomp} will be made
##' available to each of the user-defined functions.  To prevent errors due to
##' misspelling, a warning is issued if any such arguments are detected.
##' @return The \code{pomp} constructor function returns an object, call it
##' \code{P}, of class \sQuote{pomp}.  \code{P} contains, in addition to the
##' data, any elements of the model that have been specified as arguments to
##' the \code{pomp} constructor function.  One can add or modify elements of
##' \code{P} by means of further calls to \code{pomp}, using \code{P} as the
##' first argument in such calls.
##'
##' @section Important note:
##'
##' \strong{ It is not typically necessary (or even
##' feasible) to define all of the basic components for any given purpose.
##' Each \pkg{pomp} algorithm makes use of only a subset of these components.
##' Any algorithm requiring a component that is not present will generate an
##' error letting you know that you have not provided a needed component.  }
##'
##' @author Aaron A. King
##' @references
##' A. A. King, D. Nguyen, and E. L. Ionides (2016)
##' Statistical Inference for Partially Observed Markov Processes via the Package \pkg{pomp}.
##' Journal of Statistical Software 69(12): 1--43.
##'
##' D. T. Gillespie (1977)
##' Exact stochastic simulation of coupled chemical reactions.
##' Journal of Physical Chemistry 81:2340--2361.

##' @rdname pomp
pomp <- function (data, times, t0, ..., rinit, rprocess, dprocess,
  rmeasure, dmeasure, skeleton, rprior, dprior, partrans,
  params, covar, zeronames,
  obsnames, statenames, paramnames, covarnames,
  PACKAGE, globals, cdir, cfile, shlib.args) {

  ep <- paste0("in ",sQuote("pomp"),": ")

  if (missing(data))
    stop(ep,sQuote("data")," is a required argument.",call.=FALSE)

  if (!inherits(data,what=c("data.frame","pomp","NULL")))
    stop(ep,sQuote("data")," must be a data frame or an object of class ",
      sQuote("pomp"),".",call.=FALSE)

  ## return as quickly as possible if no work is to be done
  if (nargs()==1) return(data)

  construct_pomp(
    data=data,times=times,t0=t0,...,
    rinit=rinit,rprocess=rprocess,dprocess=dprocess,
    rmeasure=rmeasure,dmeasure=dmeasure,
    skeleton=skeleton,rprior=rprior,dprior=dprior,partrans=partrans,
    params=params, covar=covar,zeronames=zeronames,
    obsnames=obsnames,statenames=statenames,paramnames=paramnames,
    covarnames=covarnames,PACKAGE=PACKAGE,
    globals=globals,cdir=cdir,cfile=cfile,shlib.args=shlib.args
  )
}

setGeneric("construct_pomp",
  function(data,...)standardGeneric("construct_pomp"))

setMethod(
  "construct_pomp",
  signature=signature(data="data.frame"),
  definition = function (data, times, ...) {

    ep <- paste0("in ",sQuote("pomp"),": ")  # error prefix

    if (anyDuplicated(names(data))) {
      stop(ep,"names of data variables must be unique.", call.=FALSE)
    }
    if (missing(times))
      stop(ep,sQuote("times")," is a required argument",call.=FALSE)
    if ((is.numeric(times) && (times<1 || times>ncol(data) ||
        times!=as.integer(times))) ||
        (is.character(times) && (!(times%in%names(data)))) ||
        (!is.numeric(times) && !is.character(times)) || length(times)!=1) {
      stop(ep,"when ",sQuote("data")," is a data frame, ",sQuote("times"),
        " must identify a single column of ",sQuote("data"),
        " either by name or by index.",call.=FALSE)
    }
    if (is.numeric(times)) {
      tpos <- as.integer(times)
      timename <- names(data)[tpos]
    } else if (is.character(times)) {
      tpos <- match(times,names(data))
      timename <- times
    }
    times <- data[[tpos]]
    data <- do.call(rbind,lapply(data[-tpos],as.double))

    construct_pomp(data=data,times=times,...,timename=timename)
  }
)

setMethod(
  "construct_pomp",
  signature=signature(data="NULL"),
  definition = function (data, times, ...) {

    ep <- paste0("in ",sQuote("pomp"),": ")

    if (missing(times))
      stop(ep,sQuote("times")," is a required argument.",call.=FALSE)

    construct_pomp(data=array(dim=c(0,length(times))),times=times,...)
  }
)

setMethod(
  "construct_pomp",
  signature=signature(data="array"),
  definition = function (data, times, t0, ...,
    rinit, rprocess, dprocess, rmeasure, dmeasure, skeleton, rprior, dprior,
    partrans, params, covar, timename) {

    ep <- paste0("in ",sQuote("pomp"),": ")

    if (missing(rinit)) rinit <- NULL

    if (missing(rprocess) || is.null(rprocess)) {
      rprocess <- rproc_plugin()
    }

    if (missing(dprocess)) dprocess <- NULL
    if (missing(rmeasure)) rmeasure <- NULL
    if (missing(dmeasure)) dmeasure <- NULL

    if (missing(skeleton) || is.null(skeleton)) {
      skeleton <- skel_plugin()
    }

    if (missing(partrans) || is.null(partrans)) {
      partrans <- parameter_trans()
    }

    if (missing(rprior)) rprior <- NULL
    if (missing(dprior)) dprior <- NULL

    if (missing(params)) params <- numeric(0)
    if (is.list(params)) params <- unlist(params)

    if (missing(covar)) covar <- covariate_table()

    tryCatch(
      pomp.internal(
        data=data,
        times=times,
        t0=t0,
        timename=timename,
        rinit=rinit,
        rprocess=rprocess,
        dprocess=dprocess,
        rmeasure=rmeasure,
        dmeasure=dmeasure,
        skeleton=skeleton,
        dprior=dprior,
        rprior=rprior,
        partrans=partrans,
        params=params,
        covar=covar,
        ...
      ),
      error = function (e) {
        stop(ep,conditionMessage(e),call.=FALSE)
      }
    )
  }
)

setMethod(
  "construct_pomp",
  signature=signature(data="pomp"),
  definition = function (data, times, t0, ...,
    rinit, rprocess, dprocess, rmeasure, dmeasure, skeleton, rprior, dprior,
    partrans, params, covar, zeronames, timename) {

    ep <- paste0("in ",sQuote("pomp"),": ")  # error prefix

    if (missing(times)) {
      times <- data@times
    } else {
      time(data) <- times
    }

    if (missing(t0)) t0 <- data@t0
    if (missing(timename)) timename <- data@timename

    if (missing(rinit)) rinit <- data@rinit

    if (missing(rprocess)) {
      rprocess <- data@rprocess
    } else if (is.null(rprocess)) {
      rprocess <- rproc_plugin()
    }

    if (missing(dprocess)) dprocess <- data@dprocess
    if (missing(rmeasure)) rmeasure <- data@rmeasure
    if (missing(dmeasure)) dmeasure <- data@dmeasure

    if (missing(skeleton)) {
      skeleton <- data@skeleton
    } else if (is.null(skeleton)) {
      skeleton <- skel_plugin()
    }

    if (missing(partrans)) {
      partrans <- data@partrans
    } else if (is.null(partrans)) {
      partrans <- parameter_trans()
    }

    if (missing(rprior)) rprior <- data@rprior
    if (missing(dprior)) dprior <- data@dprior

    if (missing(params)) params <- data@params
    if (missing(covar)) covar <- data@covar
    if (missing(zeronames)) zeronames <- data@zeronames

    tryCatch(
      pomp.internal(
        data=data@data,
        times=times,
        t0=t0,
        timename=timename,
        rinit=rinit,
        rprocess=rprocess,
        dprocess=dprocess,
        rmeasure=rmeasure,
        dmeasure=dmeasure,
        skeleton=skeleton,
        rprior=rprior,
        dprior=dprior,
        partrans=partrans,
        covar=covar,
        zeronames=zeronames,
        params=params,
        .solibs=data@solibs,
        .userdata=data@userdata,
        ...
      ),
      error = function (e) {
        stop(ep,conditionMessage(e),call.=FALSE)
      }
    )
  }
)

pomp.internal <- function (data, times, t0, timename, ...,
  rinit, rprocess, dprocess, rmeasure, dmeasure, skeleton, rprior, dprior,
  partrans, params, covar, zeronames, obsnames, statenames,
  paramnames, covarnames, PACKAGE, globals, cdir, cfile, shlib.args,
  .userdata, .solibs = list(), verbose = getOption("verbose",FALSE)) {

  ep <- character(0)
  wp <- paste0("in ",sQuote("pomp"),": ")

  ## check times
  if (missing(times) || !is.numeric(times) || !all(is.finite(times)) ||
      (length(times)>1 && !all(diff(times)>0)))
    stop(ep,sQuote("times")," must be specified as an increasing sequence ",
      "of numbers.",call.=FALSE)
  storage.mode(times) <- "double"

  ## check t0
  if (missing(t0) || !is.numeric(t0) || !is.finite(t0) ||
      length(t0) > 1 ||  t0 > times[1])
    stop(ep,sQuote("t0")," must be specified as a single number not ",
      "greater than ",sQuote("times[1]"),".",call.=FALSE)
  storage.mode(t0) <- "double"

  if (missing(timename) || is.null(timename))
    timename <- "time"
  else
    timename <- as.character(timename)

  if (missing(.userdata)) .userdata <- list()
  added.userdata <- list(...)
  if (length(added.userdata)>0) {
    message(wp,"the unrecognized ",
      ngettext(length(added.userdata),"argument","arguments")," ",
      paste(sQuote(names(added.userdata)),collapse=","),
      " will be stored for use by user-defined functions."
    )
    .userdata[names(added.userdata)] <- added.userdata
  }

  if (!is(rprocess,"rprocPlugin")) {
    stop(ep,sQuote("rprocess"),
      " must be specified using one of the plugins:\n",
      sQuote("onestep.sim"),", ",sQuote("discrete.time.sim"),
      ", ",sQuote("euler.sim"),", ",sQuote("gillespie.sim"),
      ", or ",sQuote("gillespie.hl.sim"),".",call.=FALSE)
  }

  if (!is(skeleton,"skelPlugin")) {
    stop(ep,sQuote("skeleton")," must be specified using either ",
      sQuote("map")," or ",sQuote("vectorfield"),".",call.=FALSE)
  }

  if (!is(partrans,"partransPlugin")) {
    stop(ep,sQuote("partrans")," must be specified using ",
      sQuote("parameter_trans"),".",call.=FALSE)
  }

  if (missing(statenames)) statenames <- NULL
  if (missing(paramnames)) paramnames <- NULL
  if (missing(obsnames)) obsnames <- NULL
  if (missing(covarnames)) covarnames <- NULL

  statenames <- as.character(statenames)
  paramnames <- as.character(paramnames)
  obsnames <- as.character(obsnames)
  covarnames <- as.character(covarnames)

  if (missing(zeronames)) zeronames <- NULL
  zeronames <- unique(as.character(zeronames))

  ## store the data as double-precision matrix
  storage.mode(data) <- "double"
  if (length(obsnames) == 0) obsnames <- rownames(data)

  ## check the parameters and force them to be double-precision
  params <- setNames(as.double(params),names(params))
  if (length(params) > 0) {
    if (is.null(names(params)) || !is.numeric(params) ||
        !all(nzchar(names(params))))
      stop(sQuote("params")," must be a named numeric vector.",call.=FALSE)
  }

  ## use default rinit?
  default.init <- is.null(rinit) ||
    (is(rinit,"pomp_fun") && rinit@mode == pompfunmode$undef )
  if (default.init) rinit <- pomp_fun(slotname="rinit")

  if (is(rinit,"Csnippet") && length(statenames)==0) {
    stop(ep,"when ",sQuote("rinit")," is provided as a C snippet, ",
      "you must also provide ",sQuote("statenames"),".",call.=FALSE)
  }

  ## by default, use flat improper prior
  if (is.null(dprior))
    dprior <- pomp_fun(f="_pomp_default_dprior",PACKAGE="pomp")

  ## check and arrange covariates
  if (is.null(covar)) {
    covar <- covariate_table()
  } else if (!is(covar,"covartable")) {
    stop(ep,"bad option for ",sQuote("covar"),".",call.=FALSE)
  }

  if (length(covarnames) == 0) {
    covarnames <- get_covariate_names(covar)
  } else {
    covar <- tryCatch(
      select_covariates(covar,covarnames),
      error = function (e) {
        stop(ep,conditionMessage(e),call.=FALSE)
      }
    )
  }

  hitches <- hitch(
    rinit=rinit,
    step.fn=rprocess@step.fn,
    rate.fn=rprocess@rate.fn,
    dprocess=dprocess,
    rmeasure=rmeasure,
    dmeasure=dmeasure,
    rprior=rprior,
    dprior=dprior,
    toEst=partrans@to,
    fromEst=partrans@from,
    skeleton=skeleton@skel.fn,
    templates=workhorse_templates,
    obsnames=obsnames,
    statenames=statenames,
    paramnames=paramnames,
    covarnames=covarnames,
    PACKAGE=PACKAGE,
    globals=globals,
    cfile=cfile,
    cdir=cdir,
    shlib.args=shlib.args,
    verbose=verbose
  )

  ## check to make sure 'covars' is included as an argument where needed
  covar_fun_warning(covar,"rprocess",hitches$funs,wp)
  covar_fun_warning(covar,"dprocess",hitches$funs,wp)
  covar_fun_warning(covar,"rmeasure",hitches$funs,wp)
  covar_fun_warning(covar,"dmeasure",hitches$funs,wp)
  covar_fun_warning(covar,"skeleton",hitches$funs,wp)

  ## check to see if covariate times embrace the data times
  covar_time_warning(covar,times,t0,wp)

  new(
    "pomp",
    data = data,
    times = times,
    t0 = t0,
    timename = timename,
    default.init = default.init,
    rinit = hitches$funs$rinit,
    rprocess = rproc_plugin(
      rprocess,
      step.fn=hitches$funs$step.fn,
      rate.fn=hitches$funs$rate.fn
    ),
    dprocess = hitches$funs$dprocess,
    dmeasure = hitches$funs$dmeasure,
    rmeasure = hitches$funs$rmeasure,
    skeleton = skel_plugin(
      skeleton,
      skel.fn=hitches$funs$skeleton
    ),
    dprior = hitches$funs$dprior,
    rprior = hitches$funs$rprior,
    partrans = parameter_trans(
      toEst=hitches$funs$toEst,
      fromEst=hitches$funs$fromEst
    ),
    params = params,
    covar = covar,
    zeronames = zeronames,
    solibs = if (is.null(hitches$lib)) {
      .solibs
    } else {
      c(list(hitches$lib),.solibs)
    },
    userdata = .userdata
  )
}
