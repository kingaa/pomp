##' Covariates
##'
##' Including time-varying covariates in a model.
##'
##' @name covariate_table
##' @rdname covariate_table
##' @aliases covariate_table covariate_table,missing-method
##' covariate_table,ANY-method
##' @include pomp_class.R
##' @family information on model implementation
##'
##' @details
##' If the \sQuote{pomp} object contains covariates (specified via the \code{covar} argument), then interpolated values of the covariates will be available to each of the model components whenever it is called.
##' In particular, variables with names as they appear in the \code{covar} covariate table will be available to any C snippet.
##' When a basic component is defined using an \R function, that function will be called with an extra argument, \code{covars}, which will be a named numeric vector containing the interpolated values from the covariate table.
##'
##' An exception to this rule is the prior (\code{rprior} and \code{dprior}):
##' covariate-dependent priors are not allowed.
##' Nor are parameter transformations permitted to depend upon covariates.
NULL

setClass(
  "covartable",
  slots=c(
    times="numeric",
    table="matrix"
  ),
  prototype=prototype(
    times=numeric(0),
    table=array(data=numeric(0),dim=c(0,0))
  )
)

setGeneric(
  "covariate_table",
  function (..., times)
    standardGeneric("covariate_table")
)

setMethod(
  "show",
  signature=signature(object="covartable"),
  definition=function (object) {
    if (length(object@times)>0) {
      cat("  -",nrow(object@table),"records of",
        ncol(object@table),"covariates,",
        "recorded from t =",min(object@times),
        "to",max(object@times),"\n")
      cat("  - summary of covariates:\n")
      print(summary(as.data.frame(object@table)))
    } else {
      cat("\t\t<none>\n")
    }
  }
)

##' @export
setMethod(
  "covariate_table",
  signature=signature(times="missing"),
  definition=function (...) {
    ep <- paste0("in ",sQuote("covariate_table"),": ")
    if (nargs() > 0)
      stop(ep,sQuote("times")," is a required argument",call.=FALSE)
    new("covartable")
  }
)

##' @name covariate_table-ANY
##' @aliases covariate_table,ANY-method
##' @rdname covariate_table
##'
##' @param times the times corresponding to the covariates.
##' This may be given as a vector of (increasing, finite) numerical values.
##' Alternatively, one can indicate one of the variables given (either as a vector or as a data-frame column) by name or by index.
##' @param \dots numeric vectors or data frames containing time-varying covariates
##' @export
##'
setMethod(
  "covariate_table",
  signature=signature(times="ANY"),
  definition=function (..., times) {

    ep <- paste0("in ",sQuote("covariate_table"),": ")

    df <- tryCatch(
      data.frame(...,check.names=FALSE),
      error = function (e) {
        stop(ep,conditionMessage(e),call.=FALSE)
      }
    )

    if (anyDuplicated(names(df))) {
      stop(ep,"names of covariates must be unique.", call.=FALSE)
    }

    if (length(times) == 0 ||
        (length(times) > 1 && !is.numeric(times)) ||
        (length(times) == 1 && !is.numeric(times) && !is.character(times)))
      stop(ep,sQuote("times")," should either be a vector of times or ",
        "identify a single time variable.",call.=FALSE)

    if (length(times) == 1) {
      if (is.character(times)) {
        tpos <- match(times,names(df))
      } else if (is.numeric(times)) {
        tpos <- as.integer(times)
      }
      if (tpos < 1 || tpos > length(df) || !is.finite(tpos))
        stop(ep,sQuote("times")," must identify a single variable ",
          "either by name or by index.",call.=FALSE)
      times <- df[[tpos]]
      df <- df[-tpos]
    }

    if (length(df) == 0)
      stop(ep,"no covariates specified.",call.=FALSE)

    if (length(times) != nrow(df))
      stop(ep,sQuote("times")," must agree in length with the covariates.",
        call.=FALSE)

    if (any(!is.finite(times)) || !all(diff(times)>0))
      stop(ep,sQuote("times")," should be an increasing sequence of times ",
        "(without missing values).",call.=FALSE)

    new("covartable",times=as.double(times),
      table=do.call(cbind,lapply(df,as.double)))

  }
)

get_covariate_names <- function (object) {
  colnames(object@table)
}

select_covariates <- function (object, vars) {
  cnames <- colnames(object@table)
  if (!all(vars %in% cnames)) {
    missing <- vars[!(vars%in%cnames)]
    stop("variable(s) ",paste(sapply(missing,sQuote),collapse=","),
      " are not among the covariates.",call.=FALSE)
  }
  object@table <- object@table[,vars,drop=FALSE]
  object
}

covar_fun_warning <- function (object, slotname, funs, wp) {
  if (length(object@times) > 0 && !is.null(funs[[slotname]])) {
    if (funs[[slotname]]@mode == pompfunmode$Rfun &&
        !("covars" %in% names(formals(funs[[slotname]]@R.fun))))
      warning(wp,"a covariate table has been given, yet the ",
        sQuote(slotname)," function does not have ",
        sQuote("covars")," as a formal argument.",call.=FALSE)
  }
}

covar_time_warning <- function (object, times, t0, wp) {
  if ((length(object@times)>0) &&
      ((min(object@times)>t0) || (max(object@times)<max(times))))
    warning(wp,"the supplied covariate times do not embrace the ",
      "data times: covariates may be extrapolated.",call.=FALSE)
}
