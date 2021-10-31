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
    if (length(object) > 0) {
      if (!all(vapply(object,is,logical(1),"pomp"))) {
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
  contains="list",
  validity=function (object) {
    if (length(object) > 0) {
      if (!all(vapply(object,is,logical(1),"abcd_pomp"))) {
        retval <- paste0(
          "error in ",sQuote("c"),
          ": dissimilar objects cannot be combined"
        )
        return(retval)
      }
      d <- sapply(object,function(x)dim(x@traces))
      if (!all(apply(d,1,diff)==0)) {
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
  "mif2List",
  contains="list",
  validity=function (object) {
    if (length(object) > 0) {
      if (!all(vapply(object,is,logical(1),"mif2d_pomp"))) {
        retval <- paste0(
          "error in ",sQuote("c"),
          ": dissimilar objects cannot be combined"
        )
        return(retval)
      }
      d <- sapply(object,function(x)dim(x@traces))
      if (!all(apply(d,1,diff)==0)) {
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
  contains="list",
  validity=function (object) {
    if (length(object) > 0) {
      if (!all(vapply(object,is,logical(1),"pmcmcd_pomp"))) {
        retval <- paste0(
          "error in ",sQuote("c"),
          ": dissimilar objects cannot be combined"
        )
        return(retval)
      }
      d <- sapply(object,function(x)dim(x@traces))
      if (!all(apply(d,1,diff)==0)) {
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

setClass(
  "pfilterList",
  contains="list",
  validity=function (object) {
    if (length(object) > 0) {
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

setClassUnion("Pomp",c("pomp","pompList"))

setClassUnion("Pfilter",
  c("pfilterd_pomp","wpfilterd_pomp","pfilterList")
)

setClassUnion("Abc",c("abcd_pomp","abcList"))

setClassUnion("Mif2",c("mif2d_pomp","mif2List"))

setClassUnion("Pmcmc",c("pmcmcd_pomp","pmcmcList"))

setClassUnion("listie",
  members=c("pompList","abcList","mif2List","pmcmcList","pfilterList"))
