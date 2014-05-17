## this file contains short definitions of methods for the 'mif' class

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

setMethod(
          "particles",
          signature=signature(object="mif"),
          definition=function (object, Np = 1, center = coef(object),
            sd = 0, ...) {
            particles.internal(object=object,Np=Np,center=center,sd=sd,...)
          }
          )


## extract the estimated log likelihood
setMethod('logLik','mif',function(object,...)object@loglik)

## extract the convergence record
conv.rec.internal <- function (object, pars, transform = FALSE, ...) {
  if (transform) {
    pars.proper <- names(coef(object))
    pars.improper <- setdiff(colnames(object@conv.rec),pars.proper)
    retval <- cbind(
                    t(
                      partrans(
                               object,
                               params=t(object@conv.rec[,pars.proper]),
                               dir="forward"
                               )
                      ),
                    object@conv.rec[,pars.improper]
                    )
  } else {
    retval <- object@conv.rec
  }
  if (missing(pars))
    retval
  else {
    bad.pars <- setdiff(pars,colnames(retval))
    if (length(bad.pars)>0)
      stop(
           "in ",sQuote("conv.rec"),": name(s) ",
           paste(sQuote(bad.pars),collapse=","),
           " correspond to no parameter(s) in ",
           if (transform) sQuote("conv.rec(object,transform=TRUE)")
           else sQuote("conv.rec(object,transform=FALSE)"),
           call.=FALSE
           )
    retval[,pars]
  }
}

setMethod('conv.rec','mif',
          function (object, pars, transform = FALSE, ...) {
            conv.rec.internal(object=object,pars=pars,transform=transform,...)
          }
          )

## plot mif object
setMethod(
          "plot",
          "mif",
          function (x, y = NULL, ...) {
            compare.mif(x)
          }
          )

predvarplot.mif <- function (object, pars, type = 'l', mean = FALSE, ...) {
  if (!is(object,'mif'))
    stop("predvarplot error: ",sQuote("object")," must be of class ",sQuote("mif"),call.=FALSE)
  if (missing(pars))
    pars <- object@pars
  npv <- pred.var(object,pars)/(object@random.walk.sd[pars]^2)
  if (!is.null(dim(npv))) npv <- t(npv)
  if (mean && !is.null(dim(npv)))
    npv <- rowMeans(npv)
  if (!is.null(dim(npv))) {
    matplot(time(object),npv,type=type,ylab='prediction variance',xlab='time',...)
    legend(x='topright',legend=pars,col=1:length(pars),lty=1:length(pars),bty='n')
  } else {
    plot(time(object),npv,type=type,ylab='prediction variance',xlab='time',...)
  }
}
