setMethod("pfilter",
          "mif",
          function (object, Np, coef,
                    tol = 1e-17, warn = TRUE, max.fail = 0,
                    pred.mean = FALSE,
                    pred.var = FALSE,
                    filter.mean = FALSE, ...) {
            if (missing(Np))
              Np <- object@alg.pars$Np
            if (missing(coef))
              coef <- object@coef
            p <- particles(object,Np=Np,center=coef,sd=0)
            x <- init.state(object,params=p,t0=object@t0)
            pfilter(
                    as(object,'pomp'),
                    xstart=x,
                    params=p,
                    tol=tol,
                    warn=warn,
                    max.fail=max.fail,
                    pred.mean=pred.mean,
                    pred.var=pred.var,
                    filter.mean=filter.mean,
                    ...)
          }
          )
