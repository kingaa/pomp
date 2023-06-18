##' pStop, pWarn, pMess
##'
##' Custom error, warning, and message functions.
##' @name pStop
##' @rdname pStop
##' @keywords internal
##' @include package.R
##' @param who integer or character.
##' If \code{who} is an integer, it is passed to \code{\link{sys.call}} to retrieve the name of the calling function.
##' One can also pass the name of the calling function in \code{who}.
##' In either case, the name of the calling function is included in the message.
##' @param \dots message
NULL

##' @rdname pStop
pStop <- function (..., who = -1L) {
  if (is.integer(who)) {
    who <- sys.call(who)[[1]]
  }
  who <- as.character(who)
  if (length(who) > 0L)
    stop("in ",sQuote(who[1L]),": ",...,call.=FALSE)
  else
    stop(...,call.=FALSE)
}

##' @rdname pStop
pStop_ <- function (...) {
  pStop(...,who=NULL)
}

##' @rdname pStop
pWarn <- function (..., who = -1L) {
  if (is.integer(who)) {
    who <- sys.call(who)[[1]]
  }
  who <- as.character(who)
  if (length(who) > 0L)
    warning("in ",sQuote(who[1L]),": ",...,call.=FALSE)
  else
    warning(...,call.=FALSE)
}

##' @rdname pStop
pWarn_ <- function (...) {
  pWarn(...,who=NULL)
}

##' @rdname pStop
pMess <- function (..., who = -1L) {
  if (is.integer(who)) {
    who <- sys.call(who)[[1]] #nocov
  }
  who <- as.character(who)
  if (length(who) > 0L)
    message("NOTE: in ",sQuote(who[1L]),": ",...)
  else
    message("NOTE: ",...)
}

##' @rdname pStop
pMess_ <- function (...) {
  pMess(...,who=NULL)
}

undef_method <- function (method, object) {
  o <- deparse(substitute(object))
  pStop_(sQuote(method)," is undefined for ",sQuote(o)," of class ",
    sQuote(class(object)),".")
}

reqd_arg <- function (method, object) {
  if (is.null(method) || length(method)==0)
    pStop_(sQuote(object)," is a required argument.")
  else
    pStop(who=method,sQuote(object)," is a required argument.")
}

invalid_names <- function (names) {
  is.null(names) || !all(nzchar(names)) || any(is.na(names)) || anyDuplicated(names)
}
