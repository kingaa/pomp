##' Concatenate
##'
##' Methods to concatenate objects into useful listie
##'
##' @name concat
##' @rdname concat
##' @include listie.R
##' @importFrom stats setNames
##' @keywords internal
##' 
NULL

setGeneric(
  "concat",
  function (...)
    standardGeneric("concat")
)

setMethod(
  "concat",
  signature=signature(...="missing"),
  definition=function (...) {
    NULL   #nocov
  }
)

setMethod(
  "concat",
  signature=signature(...="ANY"),
  definition=function (...) {
    pStop_(sQuote("c")," is not defined for objects of mixed class.")
  }
)

flatten_list <- function (...) {
  unlist(
    lapply(
      list(...),
      \(z) {
        if (is(z,"list")) {
          setNames(as(z,"list"),names(z))
        } else {
          z
        }
      }
    )
  )
}

##' @rdname concat
setMethod(
  "concat",
  signature=signature(...="Pomp"),
  definition=function (...) {
    y <- flatten_list(...)
    setNames(
      new(
        "pompList",
        lapply(y,as,"pomp")
      ),
      names(y)
    )
  }
)

##' @rdname concat
setMethod(
  "concat",
  signature=signature(...="Pfilter"),
  definition=function (...) {
    y <- flatten_list(...)
    setNames(
      new(
        "pfilterList",
        lapply(y,as,"pfilterd_pomp")
      ),
      names(y)
    )
  }
)

##' @rdname concat
setMethod(
  "concat",
  signature=signature(...="Abc"),
  definition=function (...) {
    y <- flatten_list(...)
    setNames(
      new(
        "abcList",
        lapply(y,as,"abcd_pomp")
      ),
      names(y)
    )
  }
)

##' @rdname concat
setMethod(
  "concat",
  signature=signature(...="Mif2"),
  definition=function (...) {
    y <- flatten_list(...)
    setNames(
      new(
        "mif2List",
        lapply(y,as,"mif2d_pomp")
      ),
      names(y)
    )
  }
)

##' @rdname concat
setMethod(
  "concat",
  signature=signature(...="Pmcmc"),
  definition=function (...) {
    y <- flatten_list(...)
    setNames(
      new(
        "pmcmcList",
        lapply(y,as,"pmcmcd_pomp")
      ),
      names(y)
    )
  }
)
