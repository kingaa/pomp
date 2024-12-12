##' Saved states
##'
##' Retrieve latent state trajectories from a particle filter calculation.
##'
##' When one calls \code{\link{pfilter}} with \code{save.states="filter"} or \code{save.states="prediction"}, the latent state vector associated with each particle is saved.
##' This can be extracted by calling \code{saved_states} on the \sQuote{pfilterd.pomp} object.
##' If the filtered particles are saved, these particles are \emph{unweighted}, saved \emph{after} resampling using their normalized weights.
##' If the argument  \code{save.states="prediction"} was used, the particles correspond to simulations from \code{rprocess}, and their corresponding unnormalized weights are included in the output. 
##'
##' @name saved_states
##' @aliases saved_states,ANY-method saved_states,missing-method
##' @include pfilter.R pmcmc.R melt.R
##' @rdname saved_states
##' @family particle filter methods
##' @family extraction methods
##' @inheritParams filter_mean
##' @param format character;
##' format of the returned object (see below).
##'
##' @return According to the \code{format} argument, the saved states are returned either as a list or a data frame.
##'
##' If \code{format="data.frame"}, then the returned data frame holds the state variables and (optionally) the unnormalized log weight of each particle at each observation time.
##' The \code{.id} variable distinguishes particles.
##'
##' If \code{format="list"} and \code{\link{pfilter}} was called with \code{save.states="unweighted"} or \code{save.states="TRUE"}, the returned list contains one element per observation time.
##' Each element consists of a matrix, with one row for each state variable and one column for each particle.
##' If \code{\link{pfilter}} was called with \code{save.states="weighted"}, the list itself contains two lists:
##' the first holds the particles as above, the second holds the corresponding unnormalized log weights.
##' In particular, it has one element per observation time; each element is the vector of per-particle log weights.
##'
NULL

setGeneric(
  "saved_states",
  function (object,...) standardGeneric("saved_states")
)

setMethod(
  "saved_states",
  signature=signature(object="missing"),
  definition=function (...) {
    reqd_arg("saved_states","object")
  }
)

setMethod(
  "saved_states",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    undef_method("saved_states",object)
  }
)

##' @rdname saved_states
##' @export
setMethod(
  "saved_states",
  signature=signature(object="pfilterd_pomp"),
  definition=function (object, ...,
    format = c("list","data.frame")) {
    format <- match.arg(format)
    if (format=="list") {
      object@saved.states
    } else if ("weights" %in% names(object@saved.states)) {
      s <- melt(object@saved.states$states)
      w <- melt(object@saved.states$weights)
      s[[object@timename]] <- time(object)[as.integer(s$.L1)]
      w[[object@timename]] <- time(object)[as.integer(w$.L1)]
      w$name <- ".log.weight"
      x <- rbind(
        s[,c(object@timename,".id","name","value")],
        w[,c(object@timename,".id","name","value")]
      )
      x <- x[order(x[[object@timename]],x$.id),]
      row.names(x) <- NULL
      x
    } else {
      s <- melt(object@saved.states)
      s[[object@timename]] <- time(object)[as.integer(s$.L1)]
      s <- s[,c(object@timename,".id","name","value")]
      row.names(s) <- NULL
      s
    }
  }
)

##' @rdname saved_states
##' @export
setMethod(
  "saved_states",
  signature=signature(object="pfilterList"),
  definition=function (object, ...,
    format = c("list","data.frame")) {
    format <- match.arg(format)
    x <- lapply(object,saved_states,...,format=format)
    if (format == "data.frame") {
      x <- rbind_fill(x,.id=".L1")
    }
    x
  }
)
