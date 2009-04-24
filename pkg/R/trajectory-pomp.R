setMethod(
          "trajectory",
          "pomp",
          function (object, params, times, ...) {
            if (missing(params)) {
              params <- coef(object)
              if (length(params)==0) {
                stop("trajectory error: ",sQuote("params")," must be supplied",call.=FALSE)
              }
            }
            nrep <- NCOL(params)
            if (is.null(dim(params))) {
              params <- matrix(
                               params,
                               nrow=length(params),
                               ncol=nrep,
                               dimnames=list(
                                 names(params),
                                 NULL
                                 )
                               )
            }
            paramnames <- rownames(params)
            if (is.null(paramnames))
              stop("pfilter error: ",sQuote("params")," must have rownames",call.=FALSE)

            if (missing(times))
              times <- time(object,t0=TRUE)

            params <- as.matrix(params)
            if (missing(times))
              times <- time(object,t0=TRUE)
            x0 <- init.state(object,params=params,t0=times[1])
            x <- array(
                       dim=c(nrow(x0),nrep,length(times)),
                       dimnames=list(rownames(x0),NULL,NULL)
                       )
            switch(
                   object@skeleton.type,
                   map={                # iterate the map
                     x[,,1] <- x0
                     for (k in 2:length(times)) {
                       x[,,k] <- skeleton(
                                          object,
                                          x=x[,,k-1,drop=FALSE],
                                          t=times[k-1],
                                          params=params
                                          )
                     }
                   },
                   vectorfield={        # integrate the vectorfield
                     for (j in 1:nrep) {
                       X <- try(
                                lsoda(
                                      y=x0[,j],
                                      times=times,
                                      func=function(t,y,parms){
                                        list(
                                             skeleton(
                                                      object,
                                                      x=array(
                                                        data=y,
                                                        dim=c(length(y),1,1),
                                                        dimnames=list(names(y),NULL,NULL)
                                                        ),
                                                      t=t,
                                                      params=as.matrix(parms)
                                                      ),
                                             NULL
                                             )
                                      },
                                      parms=params[,j],
                                      ...
                                      ),
                                silent=FALSE
                                )
                       if (inherits(X,'try-error'))
                         stop("trajectory error: error in ",sQuote("lsoda"),call.=FALSE)
                       if (attr(X,'istate')[[1]]!=2)
                         warning("abnormal exit from ",sQuote("lsoda"),", istate = ",attr(X,'istate'),call.=FALSE)
                       x[,j,] <- t(X[,-1])
                     }
                   },
                   unspecified=stop("deterministic skeleton not specified")
                   )
            x
          }
          )
