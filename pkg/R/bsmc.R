## Bayesian particle filtering codes
##
## in annotation L&W AGM == Liu & West "A General Algorithm"
## 
## params = the initial particles for the parameter values
##          these are drawn from the prior distribution for the parameters
##          if a parameter is being held fixed, it is given as a row of NA's
## est = names of parameters to estimate.  Other parameters are not updated.
## Np = number of particles
## discount = delta, introduced in section 3.3 in Liu & West
## ntries = number of samplesto draw from x_t+1|x(k)_t to estimate mean of mu(k)_t+1 as in sect 2.2 Liu & West
## ntimes = number of timesteps in observation vector
## lower  = lower bounds on prior
## upper  = upper bounds on prior

setGeneric("bsmc",function(object,...)standardGeneric("bsmc"))

setMethod(
          "bsmc",
          "pomp",
          function (object, params, est,
                    smooth = 0.1,
                    ntries = 1,
                    tol = 1e-17,
                    lower = -Inf, upper = Inf,
                    seed = NULL,
                    verbose = getOption("verbose"),
                    max.fail = 0,
                    ...) {

            if (missing(seed)) seed <- NULL
            if (!is.null(seed)) {
              if (!exists(".Random.seed",where=.GlobalEnv)) { # need to initialize the RNG
                runif(1)
              }
              save.seed <- get(".Random.seed",pos=.GlobalEnv)
              set.seed(seed)
            }

            error.prefix <- paste(sQuote("bsmc"),"error: ")

            if (missing(params)) {
              if (length(object@params)>0) {
                params <- object@params
              } else {
                stop(error.prefix,sQuote("params")," must be supplied",call.=FALSE)
              }
            }

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
            prior <- params

            npars <- nrow(params)
            paramnames <- rownames(params)
            estind <- match(est,paramnames)
            npars.est <- length(estind)
            
            if (npars.est<1)
              stop(error.prefix,"no parameters to estimate",call.=FALSE)

            if (is.null(paramnames))
              stop(error.prefix,sQuote("params")," must have rownames",call.=FALSE)

            if ((length(smooth)!=1)||(smooth>1)||(smooth<=0))
              stop(error.prefix,sQuote("smooth")," must be a scalar in [0,1)",call.=FALSE)

            hsq <- smooth^2             #  see Liu & West eq(3.6) p10
            shrink <- sqrt(1-hsq)

            if (
                ((length(lower)>1)&&(length(lower)!=npars.est))||
                ((length(upper)>1)&&(length(upper)!=npars.est))
                ) {
              stop(
                   error.prefix,
                   sQuote("lower")," and ",sQuote("upper"),
                   " must each have length 1 or length equal to that of ",sQuote("est"),
                   call.=FALSE                   
                   )
            }

            for (j in seq_len(Np)) {
              if (any((params[estind,j]<lower)|(params[estind,j]>upper))) {
                ind <- which((params[estind,j]<lower)|(params[estind,j]>upper))
                stop(
                     error.prefix,
                     "parameter(s) ",paste(paramnames[estind[ind]],collapse=","),
                     " in column ",j," in ",sQuote("params"),
                     " is/are outside the box defined by ",
                     sQuote("lower")," and ",sQuote("upper"),
                     call.=FALSE
                     )
              }
            }

            xstart <- init.state(object,params=params)
            statenames <- rownames(xstart)
            nvars <- nrow(xstart)
            
            times <- time(object,t0=TRUE)
            x <- xstart

            loglik <- rep(NA,ntimes)
            eff.sample.size <- rep(NA,ntimes)
            nfail <- 0
            
            mu <- array(data=NA,dim=c(nvars,Np,1)) 
            rownames(mu) <- rownames(xstart)
            m  <- array(data=NA,dim=c(npars,Np))
            rownames(m) <- rownames(params)
            
            for (nt in seq_len(ntimes)) {
              
              ## calculate particle means ; as per L&W AGM (1)
              params.mean <- apply(params,1,mean) 
              ## calculate particle covariances : as per L&W AGM (1)
              params.var  <- cov(t(params[estind,,drop=FALSE]))

              if (verbose) {
                cat("at step",nt,"(time =",times[nt+1],")\n")
                print(
                      rbind(
                            prior.mean=params.mean[estind],
                            prior.sd=sqrt(diag(params.var))
                            )
                      )
              }

              for (j in seq_len(Np) ) {
                x.tmp <- matrix(data=x[,j],nrow=nvars,ncol=ntries)
                rownames(x.tmp) <- statenames
                p.tmp <- matrix(data=params[,j],nrow=npars,ncol=ntries)
                rownames(p.tmp) <- paramnames

                ## update mean of states at time nt as per L&W AGM (1) 
                tries <- rprocess( 
                                  object,
                                  xstart=x.tmp,
                                  times=times[c(nt,nt+1)],
                                  params=p.tmp
                                  )

                mu[,j,1] <- apply(tries[,,2,drop=FALSE],1,mean)			

                ## shrink parameters towards mean as per Liu & West eq (3.3) and L&W AGM (1)
                m[,j] <- shrink*params[,j] + (1-shrink)*params.mean 
              }
              ## evaluate probability of obervation given mean value of parameters and states (used in L&W AGM (5) below)
              g <- dmeasure( 
                            object,
                            y=object@data[,nt,drop=FALSE],
                            x=mu,
                            times=times[nt+1],
                            params=m											
                            )	
              ## sample indices -- From L&W AGM (2)
              k <- sample.int(n=Np,size=Np,replace=TRUE,prob=g)
              params <- params[,k]
              m <- m[,k]
              g <- g[k]

              ## sample new parameter vector as per L&W AGM (3) and Liu & West eq(3.2)
              pvec <- try(
                          mvtnorm::rmvnorm(
                                           n=Np,
                                           mean=rep(0,npars.est),
                                           sigma=hsq*params.var,
                                           method="eigen"
                                           ),
                          silent=FALSE
                          )
              if (inherits(pvec,"try-error"))
                stop(error.prefix,"error in ",sQuote("rmvnorm"),call.=FALSE)
              if (any(!is.finite(pvec)))
                stop(error.prefix,"extreme particle depletion",call.=FALSE)
              params[estind,] <- m[estind,]+t(pvec)

              ## sample current state vector x^(g)_(t+1) as per L&W AGM (4)
              X <- rprocess(
                            object,
                            xstart=x[,k,drop=FALSE],
                            times=times[c(nt,nt+1)],
                            params=params
                            )[,,2,drop=FALSE]

              ## evaluate likelihood of observation given X (from L&W AGM (4))
              numer <- dmeasure(
                                object,
                                y=object@data[,nt,drop=FALSE],
                                x=X,
                                times=times[nt+1],
                                params=params
                                )
              ## evaluate weights as per L&W AGM (5)
              weights <- numer/g
              
              ## apply box constraints as per the priors          
              for (j in seq_len(Np)) {
                ## the following seems problematic: will it tend to make the boundaries repellors
                if (any((params[estind,j]<lower)|(params[estind,j]>upper))) {
                  weights[j] <- 0 
                }
                ## might this rejection method be preferable?
                ## while (any((params[estind,j]<lower)|(params[estind,j]>upper))) {
                ##   ## rejection method
                ##   pvec <- try(
                ##               mvtnorm::rmvnorm(
                ##                                n=1,
                ##                                mean=rep(0,npars.est),
                ##                                sigma=hsq*params.var,
                ##                                method="eigen"
                ##                                ),
                ##               silent=FALSE
                ##               )
                ##   if (inherits(pvec,"try-error"))
                ##     stop(error.prefix,"error in ",sQuote("rmvnorm"),call.=FALSE)
                ##   if (any(!is.finite(pvec)))
                ##     stop(error.prefix,"extreme particle depletion",call.=FALSE)
                ##   params[estind,j] <- m[estind,j]+pvec[1,]
                ## }
              }

              x[,] <- X                
              
              ## test for failure to filter
              dim(weights) <- NULL
              failures <- ((weights<tol)|(!is.finite(weights))) # test for NA weights
              all.fail <- all(failures)
              if (all.fail) {                     # all particles are lost
                if (verbose) {
                  message("filtering failure at time t = ",times[nt+1])
                }
                nfail <- nfail+1
                if (nfail > max.fail)
                  stop(error.prefix,"too many filtering failures",call.=FALSE)
                loglik[nt] <- log(tol)          # worst log-likelihood
                weights <- rep(1/Np,Np)
                eff.sample.size[nt] <- 0
              } else {                  # not all particles are lost
                ## compute log-likelihood
                loglik[nt] <- log(mean(weights))
                weights[failures] <- 0
                weights <- weights/sum(weights)
                ## compute effective sample-size
                eff.sample.size[nt] <- 1/crossprod(weights)
              }

              if (verbose) {
                cat("effective sample size =",round(eff.sample.size[nt],1),"\n")
              }

              ## Matrix with samples (columns) from filtering distribution theta.t | Y.t
              if (!all.fail) {
                ## smp <- .Call(systematic_resampling,weights)
                smp <- sample.int(n=Np,size=Np,replace=TRUE,prob=weights)
                x <- x[,smp,drop=FALSE]
                params[estind,] <- params[estind,smp,drop=FALSE]
              }
              
            }
            
            if (!is.null(seed)) {
              assign(".Random.seed",save.seed,pos=.GlobalEnv)
              seed <- save.seed
            }
            
            list(
                 post=params,
                 prior=prior,
                 eff.sample.size=eff.sample.size,
                 cond.loglik=loglik,
                 smooth=smooth,
                 seed=seed,
                 nfail=nfail,
                 loglik=sum(loglik),
                 weights=weights
                 )
          }
          )
