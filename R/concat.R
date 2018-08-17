##' Concatenate
##'
##' Methods to concatenate objects into useful listies
##'
##' @name concat
##' @rdname concat
##' @include listies.R
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
    NULL
  }
)

setMethod(
  "concat",
  signature=signature(...="ANY"),
  definition=function(...) {
    stop(sQuote("c")," is not defined for objects of mixed class.",
      call.=FALSE)
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

##' @method c Pomp
##' @rdname concat
c.Pomp <- function (...) concat(...)

##' @method c Pfilter
##' @rdname concat
c.Pfilter <- function (...) concat(...)

##' @method c Abc
##' @rdname concat
c.Abc <- function (...) concat(...)

##' @method c Mif2
##' @rdname concat
c.Mif2 <- function (...) concat(...)

##' @method c Pmcmc
##' @rdname concat
c.Pmcmc <- function (...) concat(...)
