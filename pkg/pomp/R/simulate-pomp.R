## simulate a partially-observed Markov process

simulate.internal <- function (object, nsim = 1, seed = NULL, params,
                               states = FALSE, obs = FALSE,
                               times, t0, as.data.frame = FALSE, ...) {

  if (missing(times))
    times <- time(object,t0=FALSE)
  else
    times <- as.numeric(times)

  if (missing(t0))
    t0 <- timezero(object)
  else
    t0 <- as.numeric(t0)
  
  obs <- as.logical(obs)
  states <- as.logical(states)
  as.data.frame <- as.logical(as.data.frame)

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
  
  if (!obs && !states)
    object <- as(object,"pomp")

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

  if (as.data.frame) {
    if (obs && states) {
      dm <- dim(retval$obs)
      nsim <- dm[2]
      ntimes <- dm[3]
      nm <- rownames(retval$obs)
      dim(retval$obs) <- c(dm[1],prod(dm[-1]))
      rownames(retval$obs) <- nm
      dm <- dim(retval$states)
      nm <- rownames(retval$states)
      dim(retval$states) <- c(dm[1],prod(dm[-1]))
      rownames(retval$states) <- nm
      retval <- cbind(
                      as.data.frame(t(retval$obs)),
                      as.data.frame(t(retval$states))
                      )
      retval$sim <- factor(seq_len(nsim))
      retval$time <- rep(times,each=nsim)
      retval <- retval[order(retval$sim,retval$time),]
    } else if (obs || states) {
      dm <- dim(retval)
      nsim <- dm[2]
      ntimes <- dm[3]
      nm <- rownames(retval)
      dim(retval) <- c(dm[1],prod(dm[-1]))
      rownames(retval) <- nm
      retval <- as.data.frame(t(retval))
      retval$sim <- factor(seq_len(nsim))
      retval$time <- rep(times,each=nsim)
      retval <- retval[order(retval$sim,retval$time),]
    } else {
      nsim <- length(retval)
      if (nsim > 1) {
        retval <- lapply(
                         seq_len(nsim),
                         function (k) {
                           x <- as.data.frame(retval[[k]])
                           x$sim <- as.integer(k)
                           x
                         }
                         )
        retval <- do.call(rbind,retval)
        retval$sim <- factor(retval$sim)
      } else {
        retval <- as.data.frame(retval)
        retval$sim <- factor(1)
      }
    }
    
  }

  retval
}

setMethod("simulate",signature=signature(object="pomp"),definition=simulate.internal)
