## draw a set of Np particles from the user-specified distribution

particles.internal <- function (object, Np = 1, center = coef(object), sd = 0, ...) {
  if ((length(sd)==1) && (sd == 0)) {
    sd <- rep(0,length(center))
    names(sd) <- names(center)
  }
  if (is.null(names(center)) || is.null(names(sd)))
    stop("particles error: ",sQuote("center")," and ",sQuote("sd")," must have names",call.=FALSE)
  if (length(sd)!=length(center))
    stop("particles error: ",sQuote("center")," and ",sQuote("sd")," must be of equal length",call.=FALSE)
  x <- try(
           do.call(
                   object@particles,
                   c(
                     list(Np=Np,center=center,sd=sd),
                     object@userdata
                     )
                   ),
           silent=FALSE
           )
  if (inherits(x,'try-error'))
    stop("particles error: error in user-specified ",sQuote("particles")," function",call.=FALSE)
  if (
      !is.matrix(x) ||
      Np!=ncol(x) ||
      is.null(rownames(x))
      )
    stop("particles error: user ",sQuote("particles")," function must return a matrix with Np columns and rownames",call.=FALSE)
  x
}

setGeneric('particles',function(object,...)standardGeneric("particles"))

setMethod("particles",signature=signature(object="mif"),
          function (object, Np = 1, center = coef(object), sd = 0, ...) {
            particles.internal(object=object,Np=Np,center=center,sd=sd,...)
          }
          )
