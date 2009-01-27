setMethod(
          "trajectory",
          "pomp",
          function (object, params, times, ...) {
            params <- as.matrix(params)
            if (missing(times))
              times <- time(object,t0=TRUE)
            x0 <- init.state(object,params=params,t0=times[1])
            x <- array(
                       dim=c(nrow(x0),ncol(x0),length(times)),
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
                     for (j in 1:ncol(params)) {
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
