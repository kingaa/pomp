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

##' @export
setMethod(
  "show",
  signature=signature(object="unshowable"),
  definition=function (object) {
    cat("<object of class ",sQuote(as.character(class(object))),">\n",sep="")
    invisible(NULL)
  }
)

##' @export
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
  signature=signature(object="rprocPlugin"),
  definition=function (object) {
    cat("<undefined>\n")
  }
)

setMethod(
  "show",
  signature=signature(object="onestepRprocPlugin"),
  definition=function (object) {
    cat("one-step process-model simulator\n")
    cat("  - step.fun: ")
    show(object@step.fn)
  }
)

setMethod(
  "show",
  signature=signature(object="discreteRprocPlugin"),
  definition=function (object) {
    cat("discrete-time process-model simulator\n")
    cat("  - time-step =",object@delta.t,"\n")
    cat("  - step.fun: ")
    show(object@step.fn)
  }
)

setMethod(
  "show",
  signature=signature(object="eulerRprocPlugin"),
  definition=function (object) {
    cat("Euler-method process-model simulator\n")
    cat("  - time-step =",object@delta.t,"\n")
    cat("  - step.fun: ")
    show(object@step.fn)
  }
)

setMethod(
  "show",
  signature=signature(object="gillespieRprocPlugin"),
  definition=function (object) {
    cat("Gillespie-method process-model simulator\n")
    cat("  - stoichiometry matrix:\n")
    print(object@v)
    cat("  - rate.fun: ")
    show(object@rate.fn)
  }
)
