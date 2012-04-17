setGeneric("init.state",function(object,...)standardGeneric("init.state"))

## initialize the process model
setMethod(
          'init.state',
          'pomp',
          function (object, params, t0, ...) { # the package algorithms will only use these arguments
            if (missing(t0)) t0 <- object@t0
            if (missing(params)) params <- coef(object)
            .Call(do_init_state,object,params,t0)
          }
          )
