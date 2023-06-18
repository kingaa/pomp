##' Show methods
##'
##' Display the object, according to its class.
##'
##' @name show
##' @rdname show
##' @docType methods
##' @keywords internal
##' @include pomp_class.R abc.R bsmc2.R kalman.R mif2.R nlf.R pfilter.R
##' @include pmcmc.R probe.R spect.R
##' @include probe_match.R traj_match.R spect_match.R
NULL

setClassUnion("unshowable",members=c("pomp","abcd_pomp","bsmcd_pomp",
  "kalmand_pomp","mif2d_pomp","pfilterd_pomp","pmcmcd_pomp",
  "probed_pomp","spectd_pomp","probe_match_objfun","spect_match_objfun",
  "nlf_objfun","traj_match_objfun"))

##' @rdname show
##' @export
setMethod(
  "show",
  signature=signature(object="unshowable"),
  definition=function (object) {
    cat("<object of class ",sQuote(as.character(class(object))),">\n",sep="")
    invisible(NULL)
  }
)

##' @rdname show
##' @include listie.R
##' @export
setMethod(
  "show",
  signature=signature(object="listie"),
  definition=function (object) {
    y <- as(object,"list")
    names(y) <- names(object)
    show(y)
  }
)

##' @rdname show
##' @include rprocess_spec.R
##' @export
setMethod(
  "show",
  signature=signature(object="rprocPlugin"),
  definition=function (object) {
    cat("<undefined>\n")
  }
)

##' @rdname show
##' @include rprocess_spec.R
##' @export
setMethod(
  "show",
  signature=signature(object="onestepRprocPlugin"),
  definition=function (object) {
    cat("one-step process-model simulator\n")
    cat("  - step.fun: ")
    show(object@step.fn)
  }
)

##' @rdname show
##' @include rprocess_spec.R
##' @export
setMethod(
  "show",
  signature=signature(object="discreteRprocPlugin"),
  definition=function (object) {
    cat("discrete-time process-model simulator\n")
    cat("  - timestep =",object@delta.t,"\n")
    cat("  - step.fun: ")
    show(object@step.fn)
  }
)

##' @rdname show
##' @include rprocess_spec.R
##' @export
setMethod(
  "show",
  signature=signature(object="eulerRprocPlugin"),
  definition=function (object) {
    cat("Euler-method process-model simulator\n")
    cat("  - timestep =",object@delta.t,"\n")
    cat("  - step.fun: ")
    show(object@step.fn)
  }
)

##' @rdname show
##' @include rprocess_spec.R
##' @export
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

##' @rdname show
##' @include pomp_fun.R
##' @export
setMethod(
  "show",
  signature=signature("pomp_fun"),
  definition=function (object) {

    mode <- object@mode

    if (mode==pompfunmode$Rfun) { # R function

      cat("R function\n  - definition: ")
      f <- object@R.fun
      environment(f) <- globalenv()
      print(f)

    } else if (mode==pompfunmode$native) { # user supplied native code

      cat("native function\n  - name: ",sQuote(object@native.fun),"\n",sep="")
      if (length(object@PACKAGE)>0)
        cat("  - dynamically loaded from: ",sQuote(object@PACKAGE),sep="")

    } else if (mode==pompfunmode$regNative) { # built from C snippets

      cat("native function\n  - name: ",sQuote(object@native.fun),"\n",sep="")
      if (length(object@PACKAGE)>0)
        cat("  - defined by a C snippet in library ",sQuote(object@PACKAGE),sep="")

    } else {

      cat("<default>")

    }
    cat("\n")
  }
)

##' @rdname show
##' @include parameter_trans.R
##' @export
setMethod(
  "show",
  signature=signature(object="partransPlugin"),
  definition=function (object) {
    if (object@has) {
      cat("  - to estimation scale: ")
      show(object@to)
      cat("  - from estimation scale: ")
      show(object@from)
    } else {
      cat("  - to estimation scale: <identity>\n")
      cat("  - from estimation scale: <identity>\n")
    }
  }
)

##' @rdname show
##' @include covariate_table.R
##' @export
setMethod(
  "show",
  signature=signature(object="covartable"),
  definition=function (object) {
    if (length(object@times)>0) {
      cat("\n  -",ncol(object@table),"records of",
        nrow(object@table),"covariates,",
        "recorded from t =",min(object@times),
        "to",max(object@times),"\n")
      cat("  - summary of covariates:\n")
      print(summary(as.data.frame(t(object@table))))
    } else {
      cat("<none>\n")
    }
  }
)

##' @rdname show
##' @include skeleton_spec.R
##' @export
setMethod(
  "show",
  signature=signature(object="skelPlugin"),
  definition=function (object) {
    cat("<default>\n\n")
  }
)

##' @rdname show
##' @include skeleton_spec.R
##' @export
setMethod(
  "show",
  signature=signature(object="vectorfieldPlugin"),
  definition=function (object) {
    cat("vectorfield:\n  - ")
    show(object@skel.fn)
  }
)

##' @rdname show
##' @include skeleton_spec.R
##' @export
setMethod(
  "show",
  signature=signature(object="mapPlugin"),
  definition=function (object) {
    cat("map:\n")
    cat("  - timestep =",object@delta.t,"\n")
    cat("  - ")
    show(object@skel.fn)
  }
)
