## this file contains short definitions of methods for the 'abc' class

## extract the convergence record
setMethod(
          'conv.rec',
          'abc',
          function (object, pars, ...) {
            if (missing(pars)) pars <- colnames(object@conv.rec)
            object@conv.rec[,pars]
          }
          )

## plot pmcmc object
setMethod(
          "plot",
          "abc",
          function (x, y, pars, scatter = FALSE, ...) {
            if (missing(pars)) pars <- x@pars
            if (scatter) {
              pairs(conv.rec(x, pars))
            } else {
              plot.ts(conv.rec(x,pars),main="Convergence record")
            }
          }
          )

setMethod(
          "dprior",
          signature=signature(object="abc"),
          function (object, params, log = FALSE, ...) {
            do.call(
                    object@dprior,
                    list(
                         params=params,
                         hyperparams=object@hyperparams,
                         log=log
                         )
                    )
          }
          )
