setMethod(
          "pfilter",
          "mif",
          function (object, params, Np,
                    tol = 1e-17, warn = TRUE, max.fail = 0,
                    pred.mean = FALSE,
                    pred.var = FALSE,
                    filter.mean = FALSE, ...) {
            if (missing(params))
              params <- coef(object)
            if (missing(Np))
              Np <- object@alg.pars$Np
            pfilter(
                    as(object,"pomp"),
                    params=params,
                    Np=Np,
                    tol=tol,
                    warn=warn,
                    max.fail=max.fail,
                    pred.mean=pred.mean,
                    pred.var=pred.var,
                    filter.mean=filter.mean,
                    ...)
          }
          )
