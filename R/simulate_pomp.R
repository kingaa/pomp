## simulate a partially-observed Markov process

setGeneric("simulate",function(object,nsim=1,seed=NULL,...)
  standardGeneric("simulate"))

setMethod(
  "simulate",
  signature=signature(object="pomp"),
  definition=function (object, nsim = 1, seed = NULL, params,
    states = FALSE, obs = FALSE, times, t0, as.data.frame = FALSE,
    include.data = FALSE, ..., verbose = getOption("verbose", FALSE)) {

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
      ...,
      verbose=verbose
    )

  }
)

setMethod(
  "simulate",
  signature=signature(object="missing"),
  definition=function (object, nsim = 1, seed = NULL, params,
    states = FALSE, obs = FALSE, times, t0, as.data.frame = FALSE,
    include.data = FALSE, rprocess, rmeasure, obsnames, ...,
    verbose = getOption("verbose", FALSE)) {

    ep <- paste0("in ",sQuote("simulate"),": ")

    if (missing(times))
      stop(ep,sQuote("times")," is a required argument.",call.=FALSE)
    if (missing(t0))
      stop(ep,sQuote("t0")," is a required argument.",call.=FALSE)

    object <- construct_pomp(data=NULL,times=times,t0=t0,
      rprocess=rprocess,rmeasure=rmeasure,obsnames=obsnames,...,
      verbose=verbose)

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
      verbose=verbose
    )
  }
)

simulate.internal <- function (object, nsim = 1L, seed = NULL, params,
  states = FALSE, obs = FALSE, times, t0, as.data.frame = FALSE,
  include.data = FALSE, .getnativesymbolinfo = TRUE, verbose, ...) {

  ep <- paste0("in ",sQuote("simulate"),": ")

  object <- pomp(object,...,verbose=verbose)

  obs <- as.logical(obs)
  states <- as.logical(states)
  as.data.frame <- as.logical(as.data.frame)
  include.data <- as.logical(include.data)

  if (length(nsim)!=1 || !is.numeric(nsim) || !is.finite(nsim) || nsim < 1)
    stop(ep,sQuote("nsim")," must be a positive integer.",call.=FALSE)
  nsim <- as.integer(nsim)

  ## set the random seed (be very careful about this)
  seed <- as.integer(seed)
  if (length(seed)>0) {
    if (!exists(".Random.seed",envir=.GlobalEnv)) set.seed(NULL)
    save.seed <- get(".Random.seed",envir=.GlobalEnv)
    set.seed(seed)
  }

  if (missing(params)) params <- coef(object)
  if (is.list(params)) params <- unlist(params)
  if (is.null(params)) params <- numeric(0)

  params <- as.matrix(params)
  storage.mode(params) <- "double"

  if (missing(times))
    times <- time(object,t0=FALSE)
  else
    times <- as.numeric(times)

  if (missing(t0))
    t0 <- timezero(object)
  else
    t0 <- as.numeric(t0)

  if (!obs && !states)
    object <- as(object,"pomp")

  pompLoad(object,verbose=verbose)
  on.exit(pompUnload(object,verbose=verbose))

  retval <- tryCatch(
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
    error = function (e) {
      stop(ep,conditionMessage(e),call.=FALSE)
    }
  )
  .getnativesymbolinfo <- FALSE

  ## restore the RNG state
  if (length(seed) > 0) {
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
      tryCatch(
        {
          retval <- merge(od,retval,all=TRUE)
        },
        error = function (e) {
          stop(ep,"error in merging actual and simulated data.\n",
            "Check names of data, covariates, and states for conflicts.\n",
            sQuote("merge")," error message: ",conditionMessage(e),call.=FALSE)
        }
      )
    }

    retval$sim <- ordered(retval$sim)
    if (include.data) levels(retval$sim)[1L] <- "data"
    retval <- retval[order(retval$sim,retval[[object@timename]]),]

  }

  retval
}
