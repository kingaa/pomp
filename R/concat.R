##' Concatenate
##'
##' Methods to concatenate objects into useful listie
##'
##' @name concat
##' @rdname concat
##' @aliases c
##' @include listie.R
##' @importFrom stats setNames
##' @keywords internal
NULL

setGeneric(
  "concat",
  function (...)
    standardGeneric("concat")
)

setMethod(
  "concat",
  signature=signature(...="missing"),
  definition=function(...) {
    NULL   #nocov
  }
)

setMethod(
  "concat",
  signature=signature(...="ANY"),
  definition=function(...) {
    pStop_(sQuote("c")," is not defined for objects of mixed class.")
  }
)

##' @rdname concat
setMethod(
  "concat",
  signature=signature(...="Pomp"),
  definition=function (...) {
    y <- lapply(
      list(...),
      function (z) {
        if (is(z,"list"))
          setNames(as(z,"list"),names(z))
        else
          z
      }
    )
    new("pompList",unlist(y))
  }
)

##' @rdname concat
setMethod(
  "concat",
  signature=signature(...="Pfilter"),
  definition=function (...) {
    y <- lapply(
      list(...),
      function (z) {
        if (is(z,"list"))
          setNames(as(z,"list"),names(z))
        else
          z
      }
    )
    new("pfilterList",unlist(y))
  }
)

##' @rdname concat
setMethod(
  "concat",
  signature=signature(...="Abc"),
  definition=function (...) {
    y <- lapply(
      list(...),
      function (z) {
        if (is(z,"list"))
          setNames(as(z,"list"),names(z))
        else
          z
      }
    )
    new("abcList",unlist(y))
  }
)

##' @rdname concat
setMethod(
  "concat",
  signature=signature(...="Mif2"),
  definition=function (...) {
    y <- lapply(
      list(...),
      function (z) {
        if (is(z,"list"))
          setNames(as(z,"list"),names(z))
        else
          z
      }
    )
    new("mif2List",unlist(y))
  }
)

##' @rdname concat
setMethod(
  "concat",
  signature=signature(...="Pmcmc"),
  definition=function (...) {
    y <- lapply(
      list(...),
      function (z) {
        if (is(z,"list"))
          setNames(as(z,"list"),names(z))
        else
          z
      }
    )
    new("pmcmcList",unlist(y))
  }
)

##' @rdname concat
##' @export
c.Pomp <- function (...) concat(...)

##' @rdname concat
##' @export
c.Pfilter <- function (...) concat(...)

##' @rdname concat
##' @export
c.Abc <- function (...) concat(...)

##' @rdname concat
##' @export
c.Mif2 <- function (...) concat(...)

##' @rdname concat
##' @export
c.Pmcmc <- function (...) concat(...)
