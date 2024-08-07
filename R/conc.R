##' Concatenate
##'
##' @description Internal methods to concatenate objects into useful listie.
##' @details Not exported.
##' @name conc
##' @rdname conc
##' @include listie.R
##' @importFrom stats setNames
##' @keywords internal
##'
NULL

setGeneric(
  "conc",
  function (...)
    standardGeneric("conc")
)

setMethod(
  "conc",
  signature=signature(...="ANY"),
  definition=function (...) {
    if (...length()==0L) {
      cls <- "missing"
    } else {
      cls <- unique(vapply(list(...),class,character(1L)))
    }
    pStop_(
      sQuote("c"),
      " is not defined for objects of ",
      ngettext(length(cls),"class ","classes "),
      paste(sQuote(cls),collapse=", "),"."
    )
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

##' @rdname conc
setMethod(
  "conc",
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

##' @rdname conc
setMethod(
  "conc",
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

##' @rdname conc
setMethod(
  "conc",
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

##' @rdname conc
setMethod(
  "conc",
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

##' @rdname conc
setMethod(
  "conc",
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
