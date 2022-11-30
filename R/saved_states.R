##' Saved states
##'
##' Retrieve latent state trajectories from a particle filter calculation.
##'
##' When one calls \code{\link{pfilter}} with \code{save.states=TRUE}, the latent state vector associated with each particle is saved.
##' This can be extracted by calling \code{saved.states} on the \sQuote{pfilterd.pomp} object.
##' These are the \emph{unweighted} particles, saved \emph{after} resampling.
##' 
##' @name saved.states
##' @aliases saved.states,ANY-method saved.states,missing-method
##' @include pfilter.R pmcmc.R
##' @rdname saved_states
##' @family particle filter methods
##' @family extraction methods
##' @inheritParams filter.mean
##' @param format character;
##' format of the returned object (see below).
##'
##' @return According to the \code{format} argument, the saved states are returned either as a list or a data frame.
##' 
##' If \code{format="data.frame"}, then the returned data frame holds the state variables and (optionally) the log weight of each particle at each observation time.
##' The \code{.id} variable distinguishes particles.
##'
##' If \code{format="list"} and \code{\link{pfilter}} was called with \code{save.states="unweighted"} or \code{save.states="TRUE"}, the returned list contains one element per observation time.
##' Each element consists of a matrix, with one row for each state variable and one column for each particle.
##' If \code{\link{pfilter}} was called with \code{save.states="weighted"}, the list itself contains two lists:
##' the first holds the particles as above, the second holds the corresponding log weights.
##' In particular, it has one element per observation time; each element is the vector of per-particle log weights.
##'
NULL

setGeneric(
  "saved.states",
  function (object,...) standardGeneric("saved.states")
)

setMethod(
  "saved.states",
  signature=signature(object="missing"),
  definition=function (...) {
    reqd_arg("saved.states","object")
  }
)

setMethod(
  "saved.states",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    undef_method("saved.states",object)
  }
)

##' @rdname saved_states
##' @export
setMethod(
  "saved.states",
  signature=signature(object="pfilterd_pomp"),
  definition=function (object, ...,
    format = c("list","data.frame")) {
    format <- match.arg(format)
    if (format=="list") {
      object@saved.states
    } else if (length(object@saved.states)==2L) {
      s <- melt(object@saved.states$states)
      w <- melt(object@saved.states$weights)
      s$time <- time(object)[as.integer(s$L1)]
      w$time <- time(object)[as.integer(w$L1)]
      w$variable <- "weight"
      x <- rbind(
        s[,c("time",".id","variable","value")],
        w[,c("time",".id","variable","value")]
      )
      x <- x[order(x$time,x$.id),]
      rownames(x) <- NULL
      x
    } else {
      s <- melt(object@saved.states)
      s$time <- time(object)[as.integer(s$L1)]
      s <- s[,c("time",".id","variable","value")]
      rownames(s) <- NULL
      s
    }
  }
)

##' @rdname saved_states
##' @export
setMethod(
  "saved.states",
  signature=signature(object="pfilterList"),
  definition=function (object, ...) {
    lapply(object,saved.states,...)
  }
)
