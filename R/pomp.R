##' Constructor of the basic pomp object
##'
##' This function constructs a \sQuote{pomp} object, encoding a partially-observed Markov process (\acronym{POMP}) model together with a uni- or multi-variate time series.
##' As such, it is central to all the package's functionality.
##' One implements the \acronym{POMP} model by specifying some or all of its \emph{basic components}.
##' These comprise:
##' \describe{
##' \item{rinit,}{which samples from the distribution of the state process at the zero-time;}
##' \item{rprocess,}{the simulator of the unobserved Markov state process;}
##' \item{dprocess,}{the evaluator of the probability density function for transitions of the unobserved Markov state process;}
##' \item{rmeasure,}{the simulator of the observed process, conditional on the unobserved state;}
##' \item{dmeasure,}{the evaluator of the measurement model probability density function;}
##' \item{rprior,}{which samples from a prior probability distribution on the parameters;}
##' \item{dprior,}{which evaluates the prior probability density function;}
##' \item{skeleton,}{which computes the deterministic skeleton of the unobserved state process;}
##' \item{partrans,}{which performs parameter transformations.}
##' }
##' The basic structure and its rationale are described in the \emph{Journal of Statistical Software} paper, an updated version of which is to be found on the \href{https://kingaa.github.io/pomp/}{package website}.
##'
##' Each basic component is supplied via an argument of the same name.
##' These can be given in the call to \code{pomp}, or to many of the package's other functions.
##' In any case, the effect is the same: to add, remove, or modify the basic component.
##'
##' Each basic component can be furnished using C snippets, \R functions, or pre-compiled native routine available in user-provided dynamically loaded libraries.
##'
##' @name pomp
##' @rdname pomp
##' @include pomp_class.R pomp_fun.R csnippet.R safecall.R builder.R
##' @include rinit_spec.R rprocess_spec.R rmeasure_spec.R
##' @include dprocess_spec.R dmeasure_spec.R prior_spec.R
##' @include skeleton_spec.R parameter_trans.R covariate_table.R
##' @importFrom stats setNames
##'
##' @inheritParams hitch
##'
##' @param data either a data frame holding the time series data,
##' or an object of class \sQuote{pomp},
##' i.e., the output of another \pkg{pomp} calculation.
##'
##' @param times the times at which observations are made.
##' \code{times} must indicate the column of observation times by name or index.
##' The time vector must be numeric and non-decreasing.
##' Internally, \code{data} will be internally coerced to an array with storage-mode \code{double}.
##'
##' @param t0 The zero-time, i.e., the time of the initial state.
##' This must be no later than the time of the first observation, i.e., \code{t0 <= times[1]}.
##'
##' @param rinit simulator of the initial-state distribution.
##' This can be furnished either as a C snippet, an \R function, or the name of a pre-compiled native routine available in a dynamically loaded library.
##' Setting \code{rinit=NULL} sets the initial-state simulator to its default.
##' For more information, see \link[=rinit_spec]{here}.
##'
##' @param rprocess simulator of the latent state process, specified using one of the \link[=rprocess_spec]{rprocess plugins}.
##' Setting \code{rprocess=NULL} removes the latent-state simulator.
##' For more information, \link[=rprocess_spec]{see the documentation on these plugins}.
##'
##' @param dprocess optional;
##' specification of the probability density evaluation function of the unobserved state process.
##' Setting \code{dprocess=NULL} removes the latent-state density evaluator.
##' For more information, see \link[=dprocess_spec]{here}.
##'
##' @param rmeasure simulator of the measurement model, specified either as a C snippet, an \R function, or the name of a pre-compiled native routine available in a dynamically loaded library.
##' Setting \code{rmeasure=NULL} removes the measurement model simulator.
##' For more information, see \link[=rmeasure_spec]{here}.
##'
##' @param dmeasure evaluator of the measurement model density, specified either as a C snippet, an \R function, or the name of a pre-compiled native routine available in a dynamically loaded library.
##' Setting \code{dmeasure=NULL} removes the measurement density evaluator.
##' For more information, see \link[=dmeasure_spec]{here}.
##'
##' @param skeleton optional; the deterministic skeleton of the unobserved state process.
##' Depending on whether the model operates in continuous or discrete time, this is either a vectorfield or a map.
##' Accordingly, this is supplied using either the \code{\link[=skeleton_spec]{vectorfield}} or \code{\link[=skeleton_spec]{map}} fnctions.
##' For more information, see \link[=skeleton_spec]{here}.
##' Setting \code{skeleton=NULL} removes the deterministic skeleton.
##'
##' @param rprior optional; prior distribution sampler, specified either as a C snippet, an \R function, or the name of a pre-compiled native routine available in a dynamically loaded library.
##' For more information, see \link[=prior_spec]{here}.
##' Setting \code{rprior=NULL} removes the prior distribution sampler.
##'
##' @param dprior optional; prior distribution density evaluator, specified either as a C snippet, an \R function, or the name of a pre-compiled native routine available in a dynamically loaded library.
##' For more information, see \link[=prior_spec]{here}.
##' Setting \code{dprior=NULL} resets the prior distribution to its default, which is a flat improper prior.
##'
##' @param partrans optional parameter transformations, constructed using \code{\link{parameter_trans}}.
##'
##' Many algorithms for parameter estimation search an unconstrained space of parameters.
##' When working with such an algorithm and a model for which the parameters are constrained, it can be useful to transform parameters.
##' One should supply the \code{partrans} argument via a call to \code{\link{parameter_trans}}.
##' For more information, see \link[=parameter_trans]{here}.
##' Setting \code{partrans=NULL} removes the parameter transformations, i.e., sets them to the identity transformation.
##'
##' @param covar optional covariate table, constructed using \code{\link{covariate_table}}.
##'
##' If a covariate table is supplied, then the value of each of the covariates is interpolated as needed.
##' The resulting interpolated values are made available to the appropriate basic components.
##' See the documentation for \code{\link{covariate_table}} for details.
##'
##' @param params optional; named numeric vector of parameters.
##' This will be coerced internally to storage mode \code{double}.
##'
##' @param obsnames optional character vector;
##' names of the observables.
##' It is not usually necessary to specify \code{obsnames} since, by default,
##' these are read from the names of the data variables.
##'
##' @param statenames optional character vector;
##' names of the latent state variables.
##' It is typically only necessary to supply \code{statenames} when C snippets are in use.
##'
##' @param paramnames optional character vector;
##' names of model parameters.
##' It is typically only necessary to supply \code{paramnames} when C snippets are in use.
##'
##' @param covarnames optional character vector;
##' names of the covariates.
##' It is not usually necessary to specify \code{covarnames} since, by default,
##' these are read from the names of the covariates.
##'
##' @param accumvars optional character vector;
##' contains the names of accumulator variables.
##' See \link[=accumulators]{here} for a definition and discussion of accumulator variables.
##'
##' @param \dots additional arguments supply new or modify existing model characteristics or components.
##' See \code{\link{pomp}} for a full list of recognized arguments.
##'
##' When named arguments not recognized by \code{\link{pomp}} are provided, these are made available to all basic components via the so-called \dfn{userdata} facility.
##' This allows the user to pass information to the basic components outside of the usual routes of covariates (\code{covar}) and model parameters (\code{params}).
##' See the \link[=userdata]{userdata documentation here} for information on how to use this facility.
##'
##' @param verbose logical; if \code{TRUE}, diagnostic messages will be printed to the console.
##'
##' @return
##' The \code{pomp} constructor function returns an object, call it \code{P}, of class \sQuote{pomp}.
##' \code{P} contains, in addition to the data, any elements of the model that have been specified as arguments to the \code{pomp} constructor function.
##' One can add or modify elements of \code{P} by means of further calls to \code{pomp}, using \code{P} as the first argument in such calls.
##' One can pass \code{P} to most of the \pkg{pomp} package methods via their \code{data} argument.
##'
##' @section Note:
##'
##' \strong{It is not typically necessary (or indeed often feasible) to define all of the basic components for any given purpose.
##' Each \pkg{pomp} algorithm makes use of only a subset of these components.
##' Any algorithm requiring a component that is not present will generate an error letting you know that you have not provided a needed component.
##' FIXME }
##'
##' @author Aaron A. King
##'
##' @references
##' A. A. King, D. Nguyen, and E. L. Ionides (2016)
##' Statistical Inference for Partially Observed Markov Processes via the Package \pkg{pomp}.
##' Journal of Statistical Software 69(12): 1--43.
NULL

##' @rdname pomp
##' @export
pomp <- function (data, times, t0, ...,
  rinit, rprocess, dprocess, rmeasure, dmeasure,
  skeleton, rprior, dprior, partrans, covar,
  params, accumvars,
  obsnames, statenames, paramnames, covarnames,
  PACKAGE, globals, cdir, cfile, shlib.args, compile = TRUE,
  verbose = getOption("verbose", FALSE)) {

  if (missing(data))
    reqd_arg("pomp","data")

  if (!inherits(data,what=c("data.frame","pomp","NULL")))
    pStop("pomp",sQuote("data")," must be a data frame or an object of ",
      "class ",sQuote("pomp"),".")

  ## return as quickly as possible if no work is to be done
  if (is(data,"pomp") && missing(times) && missing(t0) &&
      missing(rinit) && missing(rprocess) && missing(dprocess) &&
      missing(rmeasure) && missing(dmeasure) && missing(skeleton) &&
      missing(rprior) && missing(dprior) && missing(partrans) &&
      missing(covar) && missing(params) && missing(accumvars) &&
      length(list(...)) == 0)
    return(data)

  if (missing(times)) times <- NULL

  tryCatch(
    construct_pomp(
      data=data,times=times,t0=t0,...,
      rinit=rinit,rprocess=rprocess,dprocess=dprocess,
      rmeasure=rmeasure,dmeasure=dmeasure,
      skeleton=skeleton,rprior=rprior,dprior=dprior,partrans=partrans,
      params=params,covar=covar,accumvars=accumvars,
      obsnames=obsnames,statenames=statenames,paramnames=paramnames,
      covarnames=covarnames,PACKAGE=PACKAGE,
      globals=globals,cdir=cdir,cfile=cfile,shlib.args=shlib.args,
      compile=compile,verbose=verbose
    ),
    error = function (e) pStop_(conditionMessage(e))
  )
}

setGeneric(
  "construct_pomp",
  function (data, times, ...)
    standardGeneric("construct_pomp")
)

setMethod(
  "construct_pomp",
  signature=signature(data="ANY", times="ANY"),
  definition = function (data, times, t0, ...) {
    pStop_(sQuote("times")," should either be a numeric vector of observation",
      " times or a single name identifying the column of data that represents",
      " the observation times.")
  }
)

setMethod(
  "construct_pomp",
  signature=signature(data="data.frame", times="character"),
  definition = function (data, times, t0, ...) {

    if (anyDuplicated(names(data)))
      pStop_("names of data variables must be unique.")

    if (missing(t0)) reqd_arg(NULL,"t0")

    tpos <- match(times,names(data),nomatch=0L)

    if (length(times) != 1 || tpos == 0L)
      pStop_(sQuote("times")," does not identify a single column of ",
        sQuote("data")," by name.")

    timename <- times

    times <- data[[tpos]]
    data <- do.call(rbind,lapply(data[-tpos],as.double))

    construct_pomp(data=data,times=times,t0=t0,...,timename=timename)

  }
)

setMethod(
  "construct_pomp",
  signature=signature(data="data.frame", times="numeric"),
  definition = function (data, times, t0, ...) {

    if (anyDuplicated(names(data)))
      pStop_("names of data variables must be unique.")

    if (missing(t0)) reqd_arg(NULL,"t0")

    if (length(times) != nrow(data))
      pStop_("the length of ",sQuote("times"),
        " does not match that of the data.")

    timename <- "time"
    data <- do.call(rbind,lapply(data,as.double))

    construct_pomp(data=data,times=times,t0=t0,...,timename=timename)

  }
)

setMethod(
  "construct_pomp",
  signature=signature(data="NULL", times="numeric"),
  definition = function (data, times, t0, ...) {

    if (missing(t0)) reqd_arg(NULL,"t0")

    construct_pomp(data=array(dim=c(0,length(times))),times=times,t0=t0,...)

  }
)

setMethod(
  "construct_pomp",
  signature=signature(data="array", times="numeric"),
  definition = function (data, times, ...,
    rinit, rprocess, dprocess, rmeasure, dmeasure, skeleton, rprior, dprior,
    partrans, params, covar) {

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

    if (missing(rprior)) rprior <- NULL
    if (missing(dprior)) dprior <- NULL

    if (missing(partrans) || is.null(partrans)) {
      partrans <- parameter_trans()
    }

    if (missing(params)) params <- numeric(0)
    if (is.list(params)) params <- unlist(params)

    if (missing(covar)) covar <- covariate_table()

    pomp.internal(
      data=data,
      times=times,
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
    )

  }
)

setMethod(
  "construct_pomp",
  signature=signature(data="pomp", times="numeric"),
  definition = function (data, times, ...) {
    time(data) <- times
    construct_pomp(data,times=NULL,...)
  }
)

setMethod(
  "construct_pomp",
  signature=signature(data="pomp", times="NULL"),
  definition = function (data, times, t0, timename, ...,
    rinit, rprocess, dprocess, rmeasure, dmeasure, skeleton, rprior, dprior,
    partrans, params, covar, accumvars) {

    times <- data@times
    if (missing(t0)) t0 <- data@t0
    if (missing(timename)) timename <- data@timename

    if (missing(rinit)) rinit <- data@rinit

    if (missing(rprocess)) rprocess <- data@rprocess
    else if (is.null(rprocess)) rprocess <- rproc_plugin()

    if (missing(dprocess)) dprocess <- data@dprocess
    if (missing(rmeasure)) rmeasure <- data@rmeasure
    if (missing(dmeasure)) dmeasure <- data@dmeasure

    if (missing(skeleton)) skeleton <- data@skeleton
    else if (is.null(skeleton)) skeleton <- skel_plugin()

    if (missing(rprior)) rprior <- data@rprior
    if (missing(dprior)) dprior <- data@dprior

    if (missing(partrans)) partrans <- data@partrans
    else if (is.null(partrans)) partrans <- parameter_trans()

    if (missing(params)) params <- data@params
    if (missing(covar)) covar <- data@covar
    if (missing(accumvars)) accumvars <- data@accumvars

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
      accumvars=accumvars,
      params=params,
      .solibs=data@solibs,
      .userdata=data@userdata,
      ...
    )

  }
)

pomp.internal <- function (data, times, t0, timename, ...,
  rinit, rprocess, dprocess, rmeasure, dmeasure, skeleton, rprior, dprior,
  partrans, params, covar, accumvars, obsnames, statenames,
  paramnames, covarnames, PACKAGE, globals, cdir, cfile, shlib.args,
  compile, .userdata, .solibs = list(), verbose = getOption("verbose", FALSE)) {

  ## check times
  if (missing(times) || !is.numeric(times) || !all(is.finite(times)) ||
      (length(times)>1 && !all(diff(times)>=0)))
    pStop_(sQuote("times")," must be a non-decreasing sequence of numbers.")
  storage.mode(times) <- "double"

  ## check t0
  if (missing(t0) || !is.numeric(t0) || !is.finite(t0) ||
      length(t0) > 1 ||  t0 > times[1])
    pStop_(sQuote("t0")," must be a single number not greater than ",
      sQuote("times[1]"),".")
  storage.mode(t0) <- "double"

  if (missing(timename) || is.null(timename))
    timename <- "time"
  else
    timename <- as.character(timename)

  if (missing(.userdata)) .userdata <- list()
  added.userdata <- list(...)
  if (length(added.userdata)>0) {
    message("in ",sQuote("pomp"),": the unrecognized ",
      ngettext(length(added.userdata),"argument","arguments")," ",
      paste(sQuote(names(added.userdata)),collapse=","),
      ngettext(length(added.userdata)," is"," are"),
      " available for use by the POMP basic components."
    )
    .userdata[names(added.userdata)] <- added.userdata
  }

  if (!is(rprocess,"rprocPlugin")) {
    pStop_(sQuote("rprocess"),
      " must be specified using one of the plugins:\n",
      sQuote("onestep"),", ",sQuote("discrete_time"),
      ", ",sQuote("euler"),", ",sQuote("gillespie"),
      ", or ",sQuote("gillespie_hl"),".")
  }

  if (!is(skeleton,"skelPlugin")) {
    pStop_(sQuote("skeleton")," must be specified using either ",
      sQuote("map")," or ",sQuote("vectorfield"),".")
  }

  if (!is(partrans,"partransPlugin")) {
    pStop_(sQuote("partrans")," must be specified using ",
      sQuote("parameter_trans"),".")
  }

  if (missing(statenames)) statenames <- NULL
  if (missing(paramnames)) paramnames <- NULL
  if (missing(obsnames)) obsnames <- NULL
  if (missing(covarnames)) covarnames <- NULL

  if (missing(accumvars)) accumvars <- NULL
  accumvars <- unique(as.character(accumvars))

  ## store the data as double-precision matrix
  storage.mode(data) <- "double"
  if (is.null(obsnames)) obsnames <- rownames(data)

  ## check the parameters and force them to be double-precision
  params <- setNames(as.double(params),names(params))
  if (length(params) > 0) {
    if (is.null(names(params)) || !is.numeric(params) ||
        !all(nzchar(names(params))))
      pStop_(sQuote("params")," must be a named numeric vector.")
  }

  if (is(rinit,"Csnippet") && is.null(statenames)) {
    pStop_("when ",sQuote("rinit")," is provided as a C snippet, ",
      "you must also provide ",sQuote("statenames"),".")
  }

  if (is(rmeasure,"Csnippet") && is.null(obsnames)) {
    pStop_("when ",sQuote("rmeasure")," is provided as a C snippet, ",
           "you must also provide ",sQuote("obsnames"),".")
  }

  ## check and arrange covariates
  if (is.null(covar)) {
    covar <- covariate_table()
  } else if (!is(covar,"covartable")) {
    pStop_("bad option for ",sQuote("covar"),".")
  }

  if (is.null(covarnames)) covarnames <- get_covariate_names(covar)

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
    compile=compile,
    verbose=verbose
  )

  ## check to see if covariate times embrace the data times
  covar_time_warning(covar,times,t0,"pomp")

  new(
    "pomp",
    data = data,
    times = times,
    t0 = t0,
    timename = timename,
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
    accumvars = accumvars,
    solibs = if (is.null(hitches$lib)) {
      .solibs
    } else {
      c(list(hitches$lib),.solibs)
    },
    userdata = .userdata
  )
}
