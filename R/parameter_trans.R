##' Parameter transformations
##'
##' Equipping models with parameter transformations.
##'
##' @name parameter_trans
##' @rdname parameter_trans
##' @docType methods
##' @include pomp_fun.R csnippet.R pstop.R undefined.R
##' @aliases parameter_trans parameter_trans,missing,missing-method
##' parameter_trans,ANY,ANY-method parameter_trans,ANY,missing-method
##' parameter_trans,NULL,NULL-method parameter_trans,function,function-method
##' parameter_trans,missing,ANY-method parameter_trans,pomp_fun,pomp_fun-method
##' @family implementation_info
##'
##' @param toEst,fromEst procedures that perform transformation of model parameters to and from the estimation scale, respectively.
##' These can be furnished using C snippets, \R functions, or via procedures in an external, dynamically loaded library.
##' @param log names of parameters to be log transformed.
##' @param logit names of parameters to be logit transformed.
##' @param barycentric names of parameters to be collectively transformed according to the log barycentric transformation.
##' \strong{Important note:} variables to be log-barycentrically transformed \emph{must be adjacent} in the parameter vector.
##' @param \dots ignored.
##'
##' @inheritSection pomp Note for Windows users
##' 
##' @details
##' When parameter transformations are desired, they can be integrated into the \sQuote{pomp} object via the \code{partrans} arguments using the \code{parameter_trans} function.
##' As with the other \link[=basic_components]{basic model components}, these should ordinarily be specified using C snippets.
##' When doing so, note that:
##' \enumerate{
##'   \item The parameter transformation mapping a parameter vector from the scale used by the model codes to another scale, and the inverse transformation, are specified via a call to \preformatted{parameter_trans(toEst,fromEst)}.
##'   \item The goal of these snippets is the transformation of the parameters from the natural scale to the estimation scale, and vice-versa.
##'   If \code{p} is the name of a variable on the natural scale, its value on the estimation scale is \code{T_p}.
##'   Thus the \code{toEst} snippet computes \code{T_p} given \code{p} whilst the \code{fromEst} snippet computes \code{p} given \code{T_p}.
##'   \item Time-, state-, and covariate-dependent transformations are not allowed.
##'   Therefore, neither the time, nor any state variables, nor any of the covariates will be available in the context within which a parameter transformation snippet is executed.
##' }
##'
##' These transformations can also be specified using \R functions with arguments chosen from among the parameters.
##' Such an \R function must also have the argument \sQuote{\code{...}}.
##' In this case, \code{toEst} should transform parameters from the scale that the basic components use internally to the scale used in estimation.
##' \code{fromEst} should be the inverse of \code{toEst}.
##'
##' Note that it is the user's responsibility to make sure that the transformations are mutually inverse.
##' If \code{obj} is the constructed \sQuote{pomp} object, and \code{coef(obj)} is non-empty, a simple check of this property is \preformatted{
##'   x <- coef(obj, transform = TRUE)
##'   obj1 <- obj
##'   coef(obj1, transform = TRUE) <- x
##'   identical(coef(obj), coef(obj1))
##'   identical(coef(obj1, transform=TRUE), x)}
##'
##' One can use the \code{log} and \code{logit} arguments of \code{parameter_trans} to name variables that should be log-transformed or logit-transformed, respectively.
##' The \code{barycentric} argument can name sets of parameters that should be log-barycentric transformed.
##' 
##' Note that using the \code{log}, \code{logit}, or \code{barycentric} arguments causes C snippets to be generated.
##' Therefore, you must make sure that variables named in any of these arguments are also mentioned in \code{paramnames} at the same time.
##'
##' The logit transform is defined by
##' \deqn{\mathrm{logit}(\theta)=\log\frac{\theta}{1-\theta}.}{logit(theta) = log(theta/(1-theta)).}
##'
##' The log barycentric transformation of variables \eqn{\theta_1,\dots,\theta_n}{theta1,\dots,thetan} is given by
##' \deqn{\mathrm{logbarycentric}(\theta_1,\dots,\theta_n)=\left(\log\frac{\theta_1}{\sum_i \theta_i},\dots,\log\frac{\theta_n}{\sum_i \theta_i}\right).}{logbarycentric(theta1,\dots,thetan)=(log(theta1/sum(theta)),\dots,log(thetan/sum(theta))).}
##
NULL

setClass(
  "partransPlugin",
  slots=c(
    has="logical",
    to="ANY",
    from="ANY"
  ),
  prototype=prototype(
    has=FALSE,
    to=NULL,
    from=NULL
  )
)

setMethod(
  "undefined",
  signature=signature(object="partransPlugin"),
  definition=function (object, ...) {
    undefined(object@to) || undefined(object@from)
  }
)

##' @export
setGeneric(
  "parameter_trans",
  function (toEst, fromEst, ...)
    standardGeneric("parameter_trans")
)

##' @export
setMethod(
  "parameter_trans",
  signature=signature(toEst="NULL",fromEst="NULL"),
  definition=function(toEst, fromEst, ...) {
    new("partransPlugin",has=FALSE)
  }
)

##' @export
setMethod(
  "parameter_trans",
  signature=signature(toEst="pomp_fun",fromEst="pomp_fun"),
  definition=function(toEst, fromEst, ...) {
    if (undefined(toEst) || undefined(fromEst))
      new("partransPlugin",has=FALSE)
    else
      new("partransPlugin",has=TRUE,to=toEst,from=fromEst)
  }
)

##' @name parameter_trans-Csnippet,Csnippet
##' @aliases parameter_trans,Csnippet,Csnippet-method
##' @rdname parameter_trans
##' @export
setMethod(
  "parameter_trans",
  signature=signature(toEst="Csnippet",fromEst="Csnippet"),
  definition=function(toEst, fromEst, ..., log, logit, barycentric) {
    if (missing(log) && missing(logit) && missing(barycentric))
      new("partransPlugin",has=TRUE,to=toEst,from=fromEst)
    else
      parameter_trans.internal(toEst=as(toEst,"character"),
        fromEst=as(fromEst,"character"),log=log,logit=logit,
        barycentric=barycentric)
  }
)

##' @name parameter_trans-missing,missing
##' @aliases parameter_trans,missing,missing-method
##' @rdname parameter_trans
##' @export
setMethod(
  "parameter_trans",
  signature=signature(toEst="missing",fromEst="missing"),
  definition=function(..., log, logit, barycentric) {
    if (missing(log) && missing(logit) && missing(barycentric))
      new("partransPlugin",has=FALSE)
    else
      parameter_trans.internal(toEst=NULL,fromEst=NULL,
        log=log,logit=logit,barycentric=barycentric)
  }
)

##' @name parameter_trans-character,character
##' @aliases parameter_trans,character,character-method
##' @rdname parameter_trans
##' @export
setMethod(
  "parameter_trans",
  signature=signature(toEst="character",fromEst="character"),
  definition=function(toEst, fromEst, ...) {
    new("partransPlugin",has=TRUE,to=toEst,from=fromEst)
  }
)

##' @name parameter_trans-function,function
##' @aliases parameter_trans,function,function-method
##' @rdname parameter_trans
##' @export
setMethod(
  "parameter_trans",
  signature=signature(toEst="function",fromEst="function"),
  definition=function(toEst, fromEst, ...) {
    new("partransPlugin",has=TRUE,to=toEst,from=fromEst)
  }
)

setMethod(
  "parameter_trans",
  signature=signature(toEst="ANY",fromEst="missing"),
  definition=function(toEst, fromEst, ...) {
    pStop("parameter_trans","if one of ",sQuote("toEst"),", ",
      sQuote("fromEst")," is supplied, then so must the other be.")
  }
)

setMethod(
  "parameter_trans",
  signature=signature(toEst="missing",fromEst="ANY"),
  definition=function(toEst, fromEst, ...) {
    pStop("parameter_trans","if one of ",sQuote("toEst"),", ",
      sQuote("fromEst")," is supplied, then so must the other be.")
  }
)

setMethod(
  "parameter_trans",
  signature=signature(toEst="ANY",fromEst="ANY"),
  definition=function(toEst, fromEst, ...) {
    pStop_(sQuote("parameter_trans")," not defined for arguments of class ",
      sQuote(class(toEst)),", ",sQuote(class(fromEst)),".")
  }
)

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

parameter_trans.internal <- function (toEst = NULL, fromEst = NULL,
  ..., log, logit, barycentric) {

  if (missing(log)) log <- NULL
  if (missing(logit)) logit <- NULL
  if (missing(barycentric)) barycentric <- list()

  toEst <- as.character(toEst)
  fromEst <- as.character(fromEst)

  log <- cleanForC(as.character(log))
  logit <- cleanForC(as.character(logit))
  if (is.character(barycentric)) barycentric <- list(barycentric)
  barycentric <- lapply(lapply(barycentric,as.character),cleanForC)

  out1 <- textConnection(object="fromEst",open="a",local=TRUE)
  out2 <- textConnection(object="toEst",open="a",local=TRUE)

  if (length(log) > 0) {
    tpl1 <- "\t{%v%} = exp(T_{%v%});\n"
    tpl2 <- "\tT_{%v%} = log({%v%});\n"
    for (v in log) {
      cat(file=out1,render(tpl1,v=v))
      cat(file=out2,render(tpl2,v=v))
    }
  }

  if (length(logit) > 0) {
    tpl1 <- "\t{%v%} = expit(T_{%v%});\n"
    tpl2 <- "\tT_{%v%} = logit({%v%});\n"
    for (v in logit) {
      cat(file=out1,render(tpl1,v=v))
      cat(file=out2,render(tpl2,v=v))
    }
  }

  if (length(barycentric) > 0) {
    tpl1 <- "\tfrom_log_barycentric(&{%v%},&T_{%v%},{%n%});"
    tpl2 <- "\tto_log_barycentric(&T_{%v%},&{%v%},{%n%});"
    for (b in barycentric) {
      cat(file=out1,render(tpl1,v=b[1],n=length(b)))
      cat(file=out2,render(tpl2,v=b[1],n=length(b)))
    }
  }

  close(out1)
  close(out2)
  fromEst <- paste(fromEst,collapse="\n")
  toEst <- paste(toEst,collapse="\n")

  new("partransPlugin",has=TRUE,to=Csnippet(toEst),from=Csnippet(fromEst))
}

has_partrans <- function (object) {
  object@partrans@has
}
