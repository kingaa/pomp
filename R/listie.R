##' listie
##'
##' List-like objects.
##'
##' @name listie
##' @rdname listie
##' @keywords internal
##' @include pomp_class.R
##' @include abc.R pmcmc.R mif2.R pfilter.R
NULL

setClass(
  "pompList",
  contains="list",
  validity=function (object) {
    if (length(object) > 0L) {
      if (!all(vapply(object,is,logical(1L),"pomp"))) {
        retval <- paste0(
          "error in ",sQuote("c"),
          ": dissimilar objects cannot be combined"
        )
        return(retval)
      }
    }
    TRUE
  }
)

setClass(
  "abcList",
  contains="pompList",
  validity=function (object) {
    if (length(object) > 0L) {
      if (!all(vapply(object,is,logical(1L),"abcd_pomp"))) {
        retval <- paste0(
          "error in ",sQuote("c"),
          ": dissimilar objects cannot be combined"
        )
        return(retval)
      }
      d <- sapply(object,\(x)dim(x@traces))
      if (!all(apply(d,1L,diff)==0L)) {
        retval <- paste0(
          "error in ",sQuote("c"),
          ": to be combined, ",sQuote("abcd_pomp"),
          " objects must have chains of equal length"
        )
        return(retval)
      }
    }
    TRUE
  }
)

setClass(
  "pfilterList",
  contains="pompList",
  validity=function (object) {
    if (length(object) > 0L) {
      pftypes <- vapply(object,is,logical(1L),"pfilterd_pomp")
      wftypes <- vapply(object,is,logical(1L),"wpfilterd_pomp")
      if (!all(pftypes | wftypes)) {
        retval <- paste0(
          "error in ",sQuote("c"),
          ": dissimilar objects cannot be combined"
        )
        return(retval)
      }
    }
    TRUE
  }
)

setClass(
  "mif2List",
  contains="pfilterList",
  validity=function (object) {
    if (length(object) > 0L) {
      if (!all(vapply(object,is,logical(1L),"mif2d_pomp"))) {
        retval <- paste0(
          "error in ",sQuote("c"),
          ": dissimilar objects cannot be combined"
        )
        return(retval)
      }
      d <- sapply(object,\(x)dim(x@traces))
      if (!all(apply(d,1L,diff)==0L)) {
        retval <- paste0(
          "error in ",sQuote("c"),
          ": to be combined, ",sQuote("mif2d_pomp"),
          " objects must have chains of equal length"
        )
        return(retval)
      }
    }
    TRUE
  }
)

setClass(
  "pmcmcList",
  contains="pfilterList",
  validity=function (object) {
    if (length(object) > 0L) {
      if (!all(vapply(object,is,logical(1L),"pmcmcd_pomp"))) {
        retval <- paste0(
          "error in ",sQuote("c"),
          ": dissimilar objects cannot be combined"
        )
        return(retval)
      }
      d <- sapply(object,\(x)dim(x@traces))
      if (!all(apply(d,1L,diff)==0L)) {
        retval <- paste0(
          "error in ",sQuote("c"),
          ": to be combined, ",sQuote("pmcmcd_pomp"),
          " objects must have chains of equal length"
        )
        return(retval)
      }
    }
    TRUE
  }
)

setClassUnion("Pomp",c("pomp","pompList"))

setClassUnion("Pfilter",
  c("pfilterd_pomp","wpfilterd_pomp","pfilterList")
)

setClassUnion("Abc",c("abcd_pomp","abcList"))

setClassUnion("Mif2",c("mif2d_pomp","mif2List"))

setClassUnion("Pmcmc",c("pmcmcd_pomp","pmcmcList"))

setClassUnion("listie",
  members=c("pompList","abcList","mif2List","pmcmcList","pfilterList"))
