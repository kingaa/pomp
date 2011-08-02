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
          function (object, pars, transform = FALSE, ...) {
            if (transform) {
              pars.improper <- c("loglik","nfail")
              pars.proper <- setdiff(colnames(object@conv.rec),pars.improper)
              retval <- cbind(
                              t(
                                pomp.transform(
                                               object,
                                               params=t(object@conv.rec[,pars.proper]),
                                               dir="inverse"
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
