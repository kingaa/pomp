## simulate a partially-observed Markov process

simulate.internal <- function (object, nsim = 1, seed = NULL, params,
                               states = FALSE, obs = FALSE,
                               times, t0, ...) {
  warn.condition <- (missing(t0) && ((obs || states) || (!missing(times))))
  if (warn.condition) 
    warning(
            "The default behavior of ",sQuote("simulate")," has changed.\n",
            "See the documentation (",dQuote("pomp?simulate"),") for details\n",
            "and set the ",sQuote("times")," and ",sQuote("t0"),
            " arguments appropriately to compensate.\n",
            call.=FALSE
            )

  object <- as(object,"pomp")

  if (missing(times))
    times <- time(object,t0=FALSE)
  else
    times <- as.numeric(times)

  if (missing(t0))
    t0 <- timezero(object)
  else
    t0 <- as.numeric(t0)
  
  if (missing(params))
    params <- coef(object)
  
  if (length(params)==0)
    stop("no ",sQuote("params")," specified",call.=FALSE)

  params <- as.matrix(params)

  if (!is.null(seed)) { # set the random seed (be very careful about this)
    if (!exists('.Random.seed',envir=.GlobalEnv)) runif(1)
    save.seed <- get('.Random.seed',envir=.GlobalEnv)
    set.seed(seed)
  }
  
  retval <- try(
                .Call(
                      simulation_computations,
                      object,
                      params,
                      times,
                      t0,
                      nsim,
                      obs,
                      states
                      ),
                silent=FALSE
                )

  if (inherits(retval,'try-error'))
    stop(sQuote("simulate")," error",call.=FALSE)

  if (!is.null(seed)) {                 # restore the RNG state
    assign('.Random.seed',save.seed,envir=.GlobalEnv)
  }

  retval
}

setMethod("simulate",signature=signature(object="pomp"),definition=simulate.internal)
