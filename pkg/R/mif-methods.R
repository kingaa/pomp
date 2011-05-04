## this file contains short definitions of methods for the 'mif' class

conv.rec <- function (object, ...)
  stop("function ",sQuote("conv.rec")," is undefined for objects of class ",sQuote(class(object)))
setGeneric('conv.rec')  

## extract the estimated log likelihood
setMethod('logLik','mif',function(object,...)object@loglik)

## extract the convergence record
setMethod(
          'conv.rec',
          'mif',
          function (object, pars, ...) {
            if (missing(pars)) pars <- colnames(object@conv.rec)
            object@conv.rec[,pars]
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
