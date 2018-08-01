setClassUnion("unshowable",members=c("pomp","abcd_pomp","bsmcd_pomp",
  "kalmand_pomp","mif2d_pomp","nlfd_pomp","pfilterd_pomp","pmcmcd_pomp",
  "probe_matched_pomp","probed_pomp","spect_matched_pomp","spectd_pomp",
  "traj_matched_pomp"))
setClassUnion("listies",members=c("abcList","mif2List","pmcmcList"))

setMethod(
  "print",
  signature=signature(x="unshowable"),
  definition=function (x, ...) {
    show(x)
    invisible(x)
  }
)

setMethod(
  "print",
  signature=signature(x="listies"),
  definition=function (x, ...) {
    show(x)
    invisible(x)
  }
)


setMethod(
  "show",
  signature=signature(object="unshowable"),
  definition=function (object) {
    cat("<object of class ",sQuote(as.character(class(object))),">\n",sep="")
    invisible(NULL)
  }
)

setMethod(
  "show",
  signature=signature(object="listies"),
  definition=function (object) {
    y <- as(object,"list")
    names(y) <- names(object)
    show(y)
  }
)
