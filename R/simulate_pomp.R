## simulate a partially-observed Markov process

simulate.internal <- function (object, nsim = 1, seed = NULL, params,
                               states = FALSE, obs = FALSE,
                               times, t0, as.data.frame = FALSE,
                               include.data = FALSE,
                               .getnativesymbolinfo = TRUE, ...) {
  pompLoad(object)

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
  include.data <- as.logical(include.data)

  if (missing(params))
    params <- coef(object)
  
  if (length(params)==0)
    stop("no ",sQuote("params")," specified",call.=FALSE)

  params <- as.matrix(params)

  ## set the random seed (be very careful about this)
  seed <- as.integer(seed)
  if (length(seed)>0) {
    if (!exists('.Random.seed',envir=.GlobalEnv)) set.seed(NULL)
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
                      states,
                      .getnativesymbolinfo
                      ),
                silent=FALSE
                )
  .getnativesymbolinfo <- FALSE
  
  if (inherits(retval,'try-error'))
    stop(sQuote("simulate")," error",call.=FALSE)

  ## restore the RNG state
  if (length(seed)>0) {                 
    assign('.Random.seed',save.seed,envir=.GlobalEnv)
  }

  if (as.data.frame) {
    if (obs && states) {
      dm <- dim(retval$obs)
      nsim <- dm[2L]
      nm <- rownames(retval$obs)
      dim(retval$obs) <- c(dm[1L],prod(dm[-1L]))
      rownames(retval$obs) <- nm
      dm <- dim(retval$states)
      nm <- rownames(retval$states)
      dim(retval$states) <- c(dm[1L],prod(dm[-1L]))
      rownames(retval$states) <- nm
      retval <- cbind(
                      as.data.frame(t(retval$obs)),
                      as.data.frame(t(retval$states))
                      )
      retval$sim <- seq_len(nsim)
      retval$time <- rep(times,each=nsim)
    } else if (obs || states) {
      dm <- dim(retval)
      nsim <- dm[2L]
      nm <- rownames(retval)
      dim(retval) <- c(dm[1L],prod(dm[-1L]))
      rownames(retval) <- nm
      retval <- as.data.frame(t(retval))
      retval$sim <- seq_len(nsim)
      retval$time <- rep(times,each=nsim)
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
      } else {
        retval <- as.data.frame(retval)
        retval$sim <- 1L
      }
    }

    if (include.data) {
      od <- as.data.frame(object)
      od$sim <- 0L
      retval <- merge(od,retval,all=TRUE)
    }

    retval$sim <- ordered(retval$sim)
    if (include.data) levels(retval$sim)[1L] <- "data"
    retval <- retval[order(retval$sim,retval$time),]

  }

  pompUnload(object)

  retval
}

setMethod(
          "simulate",
          signature=signature(object="pomp"),
          definition=function (object, nsim = 1, seed = NULL, params,
            states = FALSE, obs = FALSE,
            times, t0, as.data.frame = FALSE, include.data = FALSE,
            ...)
          simulate.internal(
                            object=object,
                            nsim=nsim,
                            seed=seed,
                            params=params,
                            states=states,
                            obs=obs,
                            times=times,
                            t0=t0,
                            as.data.frame=as.data.frame,
                            include.data=include.data,
                            ...
                            )
          )
