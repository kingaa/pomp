##' Create a matrix of parameters
##'
##' \code{parmat} is a utility that makes a vector of parameters suitable for
##' use in \pkg{pomp} functions.
##' @aliases parmat parmat,missing-method parmat,ANY-method
##' @author Aaron A. King
##'
setGeneric(
  "parmat",
  function (params, ...)
    standardGeneric("parmat")
)

setMethod(
  "parmat",
  signature=signature(params="missing"),
  definition=function (...) {
    reqd_arg("parmat","params")
  }
)

setMethod(
  "parmat",
  signature=signature(params="ANY"),
  definition=function (params, ...) {
    undef_method("parmat",params)
  }
)

##' @param params named numeric vector or matrix of parameters.
##' @param nrep number of replicates (columns) desired.
##' @param names optional character; column names.
##' @param ... additional arguments, currently ignored.
##' @return \code{parmat} returns a matrix consisting of \code{nrep} copies of
##' \code{params}.
##' @example examples/ricker-bifdiag.R
##' 
##' @rdname parmat
##' @export
setMethod(
  "parmat",
  signature=signature(params="numeric"),
  definition=function (params, nrep = 1, ..., names = NULL) {
    matrix(data=params,nrow=length(params),ncol=nrep,
      dimnames=list(name=names(params),.id=names))
  }
)

##' @rdname parmat
##' @export
setMethod(
  "parmat",
  signature=signature(params="array"),
  definition=function (params, nrep = 1, ..., names = NULL) {
    d <- dim(params)
    tryCatch(
      if (length(d) == 1L) {
        matrix(
          data=as.numeric(params),
          nrow=length(params),
          ncol=nrep,
          dimnames=list(
            name=names(params),
            .id=names
          )
        )
      } else if (length(d) == 2L) {
        matrix(
          data=as.numeric(params),
          nrow=nrow(params),
          ncol=ncol(params)*nrep,
          dimnames=list(
            name=rownames(params),
            .id=names
          )
        )
      } else {
        matrix(
          data=as.numeric(params),
          nrow=nrow(params),
          ncol=prod(d[-1L])*nrep,
          dimnames=list(
            name=rownames(params),
            .id=names
          )
        )
      },
      error = function (e) pStop("parmat",conditionMessage(e)),
      warning = function (e) pStop("parmat",conditionMessage(e))
    )
  }
)

##' @rdname parmat
##' @export
setMethod(
  "parmat",
  signature=signature(params="data.frame"),
  definition=function (params, nrep = 1, ...) {
    d <- dim(params)
    rn <- rownames(params)
    nrep <- as.integer(nrep)
    if (nrep > 1)
      rn <- as.character(outer(rn,seq_len(nrep),paste,sep="_"))
    rv <- array(
      data=NA_real_,
      dim=c(d[2L],nrep*d[1L]),
      dimnames=list(name=names(params),.id=rn)
    )
    tryCatch(
      for (n in names(params)) {
        if (!is.numeric(params[[n]]))
          pStop_(sQuote("params")," must contain numeric variables only.")
        rv[n,] <- rep(as.double(params[[n]]),times=nrep)
      },
      error = function (e) pStop("parmat",conditionMessage(e)),
      warning = function (e) pStop("parmat",conditionMessage(e))
    )
    rv
  }
)
