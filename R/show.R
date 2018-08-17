##' Show methods
##'
##' Display the object, according to its class.
##'
##' @name show
##' @rdname show
##' @docType methods
##' @keywords internal
##'
##' @aliases show,unshowable-method show,listies-method show,pomp_fun-method
##' show,pomp-method show,abcd_pomp-method
##' show,bsmcd_pomp-method show,kalmand_pomp-method show,mif2d_pomp-method
##' show,nlfd_pomp-method show,pfilterd_pomp-method show,pmcmcd_pomp-method
##' show,probed_pomp-method show,spectd_pomp-method
##' show,rprocPlugin-method
##' show,discreteRprocPlugin-method show,eulerRprocPlugin-method
##' show,gillespieRprocPlugin-method show,onestepRprocPlugin-method
##' show,skelPlugin-method
##' show,mapPlugin-method show,vectorfieldPlugin-method
##' show,partransPlugin-method
##' show,covartable-method
##'
##' @include pomp_class.R abc.R bsmc2.R kalman.R mif2.R nlf.R pfilter.R
##' @include pmcmc.R probe.R spect.R
##' @include probe_match.R traj_match.R spect_match.R
NULL

setClassUnion("unshowable",members=c("pomp","abcd_pomp","bsmcd_pomp",
  "kalmand_pomp","mif2d_pomp","nlfd_pomp","pfilterd_pomp","pmcmcd_pomp",
  "probed_pomp","spectd_pomp"))

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

setMethod(
  "show",
  signature=signature("pomp_fun"),
  definition=function (object) {
    mode <- object@mode
    if (mode==pompfunmode$Rfun) { # R function
      cat("\t\t")
      show(object@R.fun)
    } else if (mode==pompfunmode$native) { # user supplied native code
      cat("\t\tnative function ",sQuote(object@native.fun),sep="")
      if (length(object@PACKAGE)>0)
        cat(", dynamically loaded from ",sQuote(object@PACKAGE),sep="")
      cat("\n")
    } else if (mode==pompfunmode$regNative) { # built from C snippets
      cat("\t\tnative function ",sQuote(object@native.fun),sep="")
      if (length(object@PACKAGE)>0)
        cat(", defined by a C snippet in library ",sQuote(object@PACKAGE),sep="")
      cat("\n")
    } else {
      cat("\t\tnot specified\n")
    }
  }
)

setMethod(
  "show",
  signature=signature(object="rprocPlugin"),
  definition=function (object) {
    cat("  - undefined\n")
  }
)

setMethod(
  "show",
  signature=signature(object="onestepRprocPlugin"),
  definition=function (object) {
    cat("  - one-step process-model simulator, step.fun:\n")
    show(object@step.fn)
  }
)

setMethod(
  "show",
  signature=signature(object="discreteRprocPlugin"),
  definition=function (object) {
    cat("  - discrete-time process-model simulator, step.fun:\n")
    show(object@step.fn)
    cat("  - time-step =",object@delta.t,"\n")
  }
)

setMethod(
  "show",
  signature=signature(object="eulerRprocPlugin"),
  definition=function (object) {
    cat("  - Euler-method process-model simulator, step.fun:\n")
    show(object@step.fn)
    cat("  - time-step =",object@delta.t,"\n")
  }
)

setMethod(
  "show",
  signature=signature(object="gillespieRprocPlugin"),
  definition=function (object) {
    cat("  - Gillespie-method process-model simulator, rate.fun:\n")
    show(object@rate.fn)
    cat("  - stoichiometry matrix:\n")
    print(object@v)
    cat("\n")
  }
)

setMethod(
  "show",
  signature=signature(object="partransPlugin"),
  definition=function (object) {
    cat("  - to estimation scale:\n")
    show(object@to)
    cat("  - from estimation scale:\n")
    show(object@from)
  }
)
