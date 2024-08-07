##' Covariates
##'
##' Incorporating time-varying covariates using lookup tables.
##'
##' If the \sQuote{pomp} object contains covariates (specified via the \code{covar} argument), then interpolated values of the covariates will be available to each of the model components whenever it is called.
##' In particular, variables with names as they appear in the \code{covar} covariate table will be available to any C snippet.
##' When a basic component is defined using an \R function, that function will be called with an extra argument, \code{covars}, which will be a named numeric vector containing the interpolated values from the covariate table.
##'
##' An exception to this rule is the prior (\code{rprior} and \code{dprior}):
##' covariate-dependent priors are not allowed.
##' Nor are parameter transformations permitted to depend upon covariates.
##'
##' @name covariates
##' @rdname covariate_table
##' @aliases covariate_table covariate_table,missing-method covariate_table,ANY-method
##' @include pomp_class.R pstop.R
##' @concept covariates
##' @family implementation information
##' @family interpolation
##' @param times the times corresponding to the covariates.
##' This may be given as a vector of (non-decreasing, finite) numerical values.
##' Alternatively, one can specify by name which of the given variables is the time variable.
##' @param order the order of interpolation to be used.
##' Options are \dQuote{linear} (the default) and \dQuote{constant}.
##' Setting \code{order="linear"} treats the covariates as piecewise linear functions of time;
##' \code{order="constant"} treats them as right-continuous piecewise constant functions.
##' @param \dots numeric vectors or data frames containing time-varying covariates.
##' It must be possible to bind these into a data frame.
##'
##' @return
##' \code{covariate_table} returns a lookup table suitable for inclusion of covariates in a \sQuote{pomp} object.
##' Specifically, this is an object of class \sQuote{covartable}.
##' @section Extrapolation:
##' If \code{t} is outside the range of the lookup table, the values will be extrapolated, and a warning will be issued.
##' The type of extrapolation performed will be constant or linear according to the \code{order} flag used when creating the table.
NULL

setClass(
  "covartable",
  slots=c(
    times="numeric",
    table="matrix",
    order="integer"
  ),
  prototype=prototype(
    times=numeric(0),
    table=array(data=numeric(0),dim=c(0,0)),
    order=1L
  )
)

setGeneric(
  "covariate_table",
  function (..., times)
    standardGeneric("covariate_table")
)

setMethod(
  "covariate_table",
  signature=signature(times="missing"),
  definition=function (...) {
    if (nargs() > 0) reqd_arg("covariate_table","times")
    new("covartable")
  }
)

setMethod(
  "covariate_table",
  signature=signature(times="ANY"),
  definition=function (..., times) {
    undef_method("covariate_table",times)
  }
)

##' @rdname covariate_table
##' @export
setMethod(
  "covariate_table",
  signature=signature(times="numeric"),
  definition=function (..., order = c("linear", "constant"), times) {
    order <- match.arg(order)
    env <- parent.frame(2)
    tryCatch(
      covariate_table_internal(...,times=times,order=order,env=env),
      error = function (e) pStop(who="covariate_table",conditionMessage(e))
    )

  }
)

##' @rdname covariate_table
##' @export
setMethod(
  "covariate_table",
  signature=signature(times="character"),
  definition=function (..., order = c("linear", "constant"), times) {
    order <- match.arg(order)
    env <- parent.frame(2)
    tryCatch(
      covariate_table_internal(...,.timevar=times,order=order,env=env),
      error = function (e) pStop(who="covariate_table",conditionMessage(e))
    )

  }
)

covariate_table_internal <- function (..., times = NULL, .timevar = NULL,
  order, env) {

  d <- as.list(substitute(list(...)))[-1]
  if (length(d)==0) pStop_("no covariates specified.")
  nm <- names(d)
  if (is.null(nm)) nm <- character(length(d))
  noname <- !nzchar(nm)
  nm[noname] <- paste0(".cov.int.",seq_len(sum(noname)))
  names(d) <- nm

  e <- new.env(parent=env)

  e$times <- as.numeric(times)

  for (i in names(d)) e[[i]] <- eval(d[[i]],envir=e)

  remove(list=c("times"),pos=e)

  e <- as.list(e,all.names=TRUE)
  names(e) <- sub(".cov.int.*","",names(e))

  df <- tryCatch(
    data.frame(e,check.names=FALSE),
    error = function (e) pStop_("binding columns: ",conditionMessage(e))
  )

  if (anyDuplicated(names(df))) pStop_("names of covariates must be unique.")

  if (!is.null(times) && (length(times) != nrow(df)))
    pStop_(sQuote("times")," must agree in length with the covariates.")

  if (!is.null(.timevar)) {

    tpos <- match(.timevar,names(df),nomatch=0L)
    if (tpos == 0L)
      pStop_(sQuote("times")," does not identify a unique time variable.")
    times <- df[[tpos]]
    df <- df[-tpos]

  }

  if (length(df) == 0) pStop_("no covariates specified.")

  if (any(!is.finite(times)) || !all(diff(times)>=0))
    pStop_(sQuote("times"),
      " must be a non-decreasing numeric sequence (without missing values).")

  names(df) <- cleanForC(names(df))

  new("covartable",times=as.double(times),
    table=do.call(rbind,lapply(df,as.double)),
    order=switch(order,linear=1L,constant=0L))

}

get_covariate_names <- function (object) {
  rownames(object@table)
}

covar_time_warning <- function (object, times, t0, wp) {
  if ((length(object@times)>0) &&
        ((min(object@times)>t0) || (max(object@times)<max(times))))
    pWarn(who=wp,"the supplied covariate times do not embrace the ",
      "data times: covariates may be extrapolated.")
}

setMethod(
  "undefined",
  signature=signature(object="covartable"),
  definition=function (object) {
    nrow(object@table) == 0L
  }
)
