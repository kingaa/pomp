## simulate a partially-observed Markov process

setMethod(
          "simulate",
          "pomp",
          function (object, nsim = 1, seed = NULL, params,
                    states = FALSE, obs = FALSE,
                    times = c(object@t0,time(object)), ...) {
            ntimes <- length(times)
            times <- as.numeric(times)
            if (ntimes<1)
              stop("if length of 'times' is less than 1, there is no work to do",call.=FALSE)
            if (missing(params)) {
              if (length(object@params)>0)
                params <- object@params
              else
                stop("no 'params' specified",call.=FALSE)
            }
            if (is.null(dim(params)))
              params <- matrix(params,ncol=1,dimnames=list(names(params),NULL))
            if (is.null(rownames(params)))
              stop("'params' must have rownames",call.=FALSE)
            npars <- ncol(params)
            if (!is.null(seed)) { # set the random seed (be very careful about this)
              if (!exists('.Random.seed',envir=.GlobalEnv)) runif(1)
              save.seed <- get('.Random.seed',envir=.GlobalEnv)
              set.seed(seed)
            }
            xstart <- init.state(object,params=params,t0=times[1])
            nreps <- npars*nsim                   # total number of replicates
            ## we will do the process model simulations with single calls to the user functions
            if (nsim > 1) {    # make nsim copies of the IC and parameter matrices
              xstart <- array(rep(xstart,nsim),dim=c(nrow(xstart),nreps),dimnames=list(rownames(xstart),NULL))
              params <- array(rep(params,nsim),dim=c(nrow(params),nreps),dimnames=list(rownames(params),NULL))
            }
            x <- rprocess(object,xstart=xstart,times=times,params=params) # simulate the process model
            ## rprocess returns a rank-3 matrix (nvars x nreps x ntimes)  
            if (obs || !states) { # we need only simulate the observation process if obs=T or states=F
              y <- rmeasure(object,x=x,times=times,params=params)
              ## rmeasure returns a rank-3 matrix (nobs x nreps x ntimes)
            }
            if (!is.null(seed)) {                 # restore the RNG state
              assign('.Random.seed',save.seed,envir=.GlobalEnv)
            }
            if (!obs && !states) { # both obs=F and states=F, return a list of pomp objects
              nm.y <- rownames(y)
              nobs <- nrow(y)
              retval <- lapply(
                               seq(length=nreps),
                               function (rep) {
                                 po <- as(object,'pomp') # copy the pomp object
                                 po@t0 <- times[1]
                                 po@times <- times[-1]
                                 po@data <- array(0,dim=c(nobs,ntimes-1),dimnames=list(nm.y,NULL))
                                 po@data[,] <- y[,rep,-1] # replace the data
                                 po@params <- params[,rep] # store the parameters
                                 po@states <- x[,rep,] # store the states
                                 po
                               }
                               )
            } else {
              if (!obs && states) {               # return states only
                retval <- x
              } else if (obs && !states) {        # return only observations
                retval <- y
              } else {                     # return both states and observations
                retval <- list(
                               states = x,
                               obs = y
                               )
              }
            }
            retval
          }
          )
