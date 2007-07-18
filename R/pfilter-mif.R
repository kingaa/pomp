pfilter.mif <- function (object,
                         Np = object@alg.pars$Np,
                         coef = object@coef,
                         ...
                         ) {
  x <- particles(object,Np=Np,center=coef,sd=0)
  pfilter(as(object,'pomp'),xstart=x$states,params=x$params,...)
}

setMethod('pfilter','mif',pfilter.mif)
