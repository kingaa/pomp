setClassUnion("unshowable",members=c("pomp","abc","bsmcd.pomp",
  "kalmand.pomp","mif2d.pomp","nlfd.pomp","pfilterd.pomp","pmcmc",
  "probe.matched.pomp","probed.pomp","spect.matched.pomp","spect.pomp",
  "traj.matched.pomp"))

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
