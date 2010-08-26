setMethod(
          "pfilter",
          "mif",
          function (object, params, Np,
                    tol = 1e-17, max.fail = 0,
                    pred.mean = FALSE,
                    pred.var = FALSE,
                    filter.mean = FALSE,
                    ...) {
            if (missing(params))
              params <- coef(object)
            if (missing(Np))
              Np <- object@alg.pars$Np
            if ("save.states"%in%names(list(...)))
              stop(
                   "pfilter error: you cannot set ",sQuote("save.states"),
                   " when ",sQuote("object")," is of class mif",
                   call.=FALSE
                   )
            pfilter(
                    object=as(object,"pomp"),
                    params=params,
                    Np=Np,
                    tol=tol,
                    max.fail=max.fail,
                    pred.mean=pred.mean,
                    pred.var=pred.var,
                    filter.mean=filter.mean,
                    save.states=FALSE,
                    ...
                    )
          }
          )
