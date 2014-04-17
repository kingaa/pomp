## this file contains short definitions of methods for the 'mif2d.pomp' class

setMethod('logLik','mif2d.pomp',function(object,...)object@loglik)

setMethod('conv.rec','mif2d.pomp',
          function (object, pars, transform = FALSE, ...) {
            conv.rec.internal(object=object,pars=pars,transform=transform,...)
          }
          )

setMethod(
          "plot",
          "mif2d.pomp",
          function (x, y = NULL, ...) {
            compare.mif(x)
          }
          )
