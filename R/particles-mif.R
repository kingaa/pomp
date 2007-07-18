## draw a set of Np particles from the user-specified distribution
particles.mif <- function (object, Np = 1, center = coef(object), sd = 0, ...) {
  if ((length(sd)==1) && (sd == 0)) {
    sd <- rep(0,length(center))
    names(sd) <- names(center)
  }
  if (is.null(names(center)) || is.null(names(sd)))
    stop("particles error: 'center' and 'sd' must have names")
  if (length(sd)!=length(center))
    stop("particles error: 'center' and 'sd' must be of equal length")
  x <- try(
           do.call(
                   object@particles,
                   c(
                     list(Np=Np,center=center,sd=sd),
                     object@userdata
                     )
                   ),
           silent=T
           )
  if (inherits(x,'try-error'))
    stop("particles error: error in user-specified 'particles' function\n",x)
  if (
      !is.list(x) ||
      !all(c('states','params')%in%names(x)) ||
      !is.matrix(x$states) ||
      !is.matrix(x$params) ||
      Np!=ncol(x$states) ||
      Np!=ncol(x$params)
      )
    stop("the 'particles' function must return a list with 'states' and 'params' matrices, each of Np columns")
  x
}


setMethod('particles','mif',particles.mif)
