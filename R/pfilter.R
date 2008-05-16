## particle filtering codes

setMethod(
          "pfilter",
          "pomp",
          function (object, params, Np,
                    tol = 1e-17, warn = TRUE, max.fail = 0,
                    pred.mean = FALSE,
                    pred.var = FALSE,
                    filter.mean = FALSE,
                    .rw.sd, ...) {
            if (missing(params)) {
              if (length(object@params)>0) {
                params <- object@params
              } else {
                stop("pfilter error: 'params' must be supplied",call.=FALSE)
              }
            }
            if (missing(Np))
              Np <- NCOL(params)
            ntimes <- length(time(object))
            if (is.null(dim(params))) {
              params <- matrix(
                               params,
                               nrow=length(params),
                               ncol=Np,
                               dimnames=list(
                                 names(params),
                                 NULL
                                 )
                               )
            }
            npars <- nrow(params)
            paramnames <- rownames(params)
            if (is.null(paramnames))
              stop("pfilter error: 'params' must have rownames",call.=FALSE)

            xstart <- init.state(object,params=params)
            statenames <- rownames(xstart)
            nvars <- nrow(xstart)

            random.walk <- !missing(.rw.sd)
            if (random.walk) {
              rw.names <- names(.rw.sd)
              if (is.null(rw.names)||!is.numeric(.rw.sd))
                stop("pfilter error: '.rw.sd' must be a named vector",call.=FALSE)
              if (any(!(rw.names%in%paramnames)))
                stop("pfilter error: the rownames of 'params' must include all of the names of '.rw.sd'",call.=FALSE)
              sigma <- .rw.sd
            }
            
            loglik <- rep(NA,ntimes)
            eff.sample.size <- rep(NA,ntimes)
            nfail <- 0
            
            pred.m <-  if (pred.mean)
              matrix(
                     data=0,
                     nrow=nvars+npars,
                     ncol=ntimes,
                     dimnames=list(c(statenames,paramnames),NULL)
                     )
            else NULL
            
            pred.v <- if (pred.var)
              matrix(
                     data=0,
                     nrow=nvars+npars,
                     ncol=ntimes,
                     dimnames=list(c(statenames,paramnames),NULL)
                     )
            else NULL
            
            filt.m <- if (filter.mean)
              matrix(
                     data=0,
                     nrow=nvars+npars,
                     ncol=ntimes,
                     dimnames=list(c(statenames,paramnames),NULL)
                     )
            else NULL

            times <- time(object,t0=TRUE)
            x <- xstart
            
            for (nt in seq(length=ntimes)) {
              
              ## advance the state variables according to the process model
              X <- try(
                       rprocess(
                                object,
                                x=x,
                                times=times[c(nt,nt+1)],
                                params=params
                                )[,,2,drop=FALSE],
                       silent=FALSE
                       )
              if (inherits(X,'try-error'))
                stop("pfilter error: process simulation error",call.=FALSE)

              x[,] <- X                 # ditch the third dimension
              
              ## prediction means
              if (pred.mean) {                    
                xx <- try(
                          c(
                            apply(x,1,mean),
                            apply(params,1,mean)
                            ),
                          silent=FALSE
                          )
                if (inherits(xx,'try-error')) {
                  stop("pfilter error: error in prediction mean computation",call.=FALSE)
                } else {
                  pred.m[,nt] <- xx
                }
              }

              ## prediction variances
              if (pred.var) {
                problem.indices <- unique(which(!is.finite(x),arr.ind=TRUE)[,1])
                if (length(problem.indices)>0) {
                  stop(
                       "pfilter error: non-finite state variables ",
                       paste(rownames(x)[problem.indices],collapse=','),
                       call.=FALSE
                       )
                }
                problem.indices <- unique(which(!is.finite(params),arr.ind=TRUE)[,1])
                if (length(problem.indices)>0) {
                  stop(
                       "pfilter error: non-finite parameters ",
                       paste(rownames(params)[problem.indices],collapse=','),
                       call.=FALSE
                       )
                }
                xx <- try(
                          c(
                            apply(x,1,var),
                            apply(params,1,var)
                            ),
                          silent=FALSE
                          )
                if (inherits(xx,'try-error')) {
                  stop("pfilter error: error in prediction variance computation",call.=FALSE)
                } else {
                  pred.v[,nt] <- xx
                }
              }

              ## determine the weights
              weights <- try(
                             dmeasure(
                                      object,
                                      y=object@data[,nt,drop=FALSE],
                                      x=X,
                                      times=times[nt+1],
                                      params=params
                                      ),
                             silent=FALSE
                             )
              if (inherits(weights,'try-error'))
                stop("pfilter error: error in calculation of weights",call.=FALSE)
              if (any(is.na(weights))) {
                ## problem.indices <- which(is.na(weights))
                stop("pfilter error: dmeasure returns NA",call.=FALSE)
              }

              ## test for failure to filter
              dim(weights) <- NULL
              failures <- weights < tol
              all.fail <- all(failures)
              if (all.fail) {                     # all particles are lost
                if (warn)
                  message("filtering failure at time t = ",times[nt+1])
                nfail <- nfail+1
                if (nfail > max.fail)
                  stop('pfilter error: too many filtering failures',call.=FALSE)
                loglik[nt] <- log(tol)          # worst log-likelihood
                weights <- rep(1/Np,Np)
                eff.sample.size[nt] <- 0
              } else {                  # not all particles are lost
                ## compute log-likelihood
                loglik[nt] <- log(mean(weights))  
                weights[failures] <- 0
                weights <- weights/sum(weights)
                ## compute effective sample-size
                eff.sample.size[nt] <- 1/(weights%*%weights) 
              }

              ## compute filtering means
              if (filter.mean) {
                filt.m[statenames,nt] <- x %*% weights
                filt.m[paramnames,nt] <- params %*% weights
              }

              ## Matrix with samples (columns) from filtering distribution theta.t | Y.t
              if (!all.fail) {
                sample <- .Call(systematic_resampling,weights)
                x <- x[,sample,drop=FALSE]
                params <- params[,sample,drop=FALSE]
              }

              ## random walk for parameters
              if (random.walk) {
                pred.v[rw.names,nt] <- pred.v[rw.names,nt]+sigma^2
                params[rw.names,] <- params[rw.names,]+rnorm(n=Np*length(sigma),mean=0,sd=sigma)
              }
            }

            list(
                 pred.mean=pred.m,
                 pred.var=pred.v,
                 filter.mean=filt.m,
                 eff.sample.size=eff.sample.size,
                 cond.loglik=loglik,
                 nfail=nfail,
                 loglik=sum(loglik)
                 )
          }
          )
