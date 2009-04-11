mif.cooling <- function (factor, n) {
  alpha <- factor^(n-1)
  list(alpha=alpha,gamma=alpha^2)
}

powerlaw.cooling <- function (init = 1, delta = 0.1, eps = (1-delta)/2, n) {
  m <- init
  if (n <= m) {                         # linear cooling regime
    alpha <- (m-n+1)/m
    gamma <- alpha
  } else {                              # power-law cooling regime
    alpha <- ((n/m)^(delta+eps))/n
    gamma <- ((n/m)^(delta+1))/n/n
  }
  list(alpha=alpha,gamma=gamma)
}

setMethod(
          "mif",
          "pomp",
          function (object, Nmif = 1,
                    start,
                    pars, ivps = character(0),
                    particles,
                    rw.sd, alg.pars,
                    weighted = TRUE, tol = 1e-17, warn = TRUE, max.fail = 0,
                    verbose = FALSE)
          {
            if (missing(rw.sd))
              stop("mif error: ",sQuote("rw.sd")," must be specified",call.=FALSE)
            if (missing(alg.pars))
              stop("mif error: ",sQuote("alg.pars")," must be specified",call.=FALSE)
            if (missing(particles)) {         # use default: normal distribution
              particles <- function (Np, center, sd, ...) {
                matrix(
                       data=rnorm(
                         n=Np*length(center),
                         mean=center,
                         sd=sd
                         ),
                       nrow=length(center),
                       ncol=Np,
                       dimnames=list(
                         names(center),
                         NULL
                         )
                       )
              }
            } else {
              particles <- match.fun(particles)
            }
            if (!all(c('Np','center','sd','...')%in%names(formals(particles))))
              stop(
                   "mif error: ",
                   sQuote("particles"),
                   " must be a function of prototype ",
                   sQuote("particles(Np,center,sd,...)"),
                   call.=FALSE
                   )
            if (missing(start)) {
              start <- coef(object)
              if (length(start)==0)
                stop("mif error: ",sQuote("start")," must be specified",call.=FALSE)
            }
            start.names <- names(start)
            if (is.null(start.names))
              stop("mif error: ",sQuote("start")," must be a named vector",call.=FALSE)

            rw.names <- names(rw.sd)
            if (is.null(rw.names) || any(rw.sd<0))
              stop("mif error: ",sQuote("rw.sd")," must be a named non-negative numerical vector",call.=FALSE)
            if (!all(rw.names%in%start.names))
              stop("mif error: all the names of ",sQuote("rw.sd")," must be names of ",sQuote("start"),call.=FALSE)

            rw.names <- names(rw.sd[rw.sd>0])
            if (missing(pars)) {
              pars <- rw.names[!(rw.names%in%ivps)]
            }
            if (length(pars) == 0)
              stop("mif error: ",sQuote("pars")," must be a nonempty character vector",call.=FALSE)
            if (
                !is.character(pars) ||
                !is.character(ivps) ||
                !all(pars%in%start.names) ||
                !all(ivps%in%start.names) ||
                any(pars%in%ivps) ||
                any(ivps%in%pars) ||
                !all(pars%in%rw.names) ||
                !all(ivps%in%rw.names)
                )
              stop(
                   "mif error: ",
                   sQuote("pars")," and ",sQuote("ivps"),
                   " must be mutually disjoint subsets of ",
                   sQuote("names(start)"),
                   " and must have a positive random-walk SDs specified in ",
                   sQuote("rw.sd"),
                   call.=FALSE
                   )

            if (!all(rw.names%in%c(pars,ivps))) {
              extra.rws <- rw.names[!(rw.names%in%c(pars,ivps))]
              warning(
                      "mif warning: the variable(s) ",
                      paste(extra.rws,collapse=", "),
                      " have positive random-walk SDs specified, but are included in neither ",
                      sQuote("pars")," nor ",sQuote("ivps"),
                      ". These random walk SDs are ignored.",
                      call.=FALSE
                      )
            }
            rw.sd <- rw.sd[c(pars,ivps)]
            rw.names <- names(rw.sd)
            
            if (!all(c('Np','cooling.factor','ic.lag','var.factor')%in%names(alg.pars)))
              stop(
                   "mif error: ",sQuote("alg.pars"),
                   " must be a named list with elements ",
                   sQuote("Np"),",",sQuote("cooling.factor"),",",sQuote("ic.lag"),
                   ",and ",sQuote("var.factor"),
                   call.=FALSE
                   )

            coef(object) <- start

            newmif <- new(
                          "mif",
                          object,
                          ivps=ivps,
                          pars=pars,
                          Nmif=as.integer(0),
                          particles=particles,
                          alg.pars=alg.pars,
                          random.walk.sd=rw.sd,
                          pred.mean=matrix(NA,0,0),
                          pred.var=matrix(NA,0,0),
                          filter.mean=matrix(NA,0,0),
                          conv.rec=matrix(
                            c(NA,NA,start),
                            nrow=1,
                            ncol=2+length(start),
                            dimnames=list(
                              0,
                              c('loglik','nfail',start.names)
                              )
                            ),
                          cond.loglik=numeric(0),
                          eff.sample.size=numeric(0),
                          loglik=as.numeric(NA)
                          )

            if (Nmif > 0) {
              mif(
                  newmif,
                  Nmif=Nmif,
                  weighted=weighted,
                  tol=tol,
                  warn=warn,
                  max.fail=max.fail,
                  verbose=verbose,
                  .ndone=0
                  )
            } else {
              newmif
            }
          }
          )
          

setMethod(
          "mif",
          "mif",
          function (object, Nmif, start, pars, ivps, rw.sd, alg.pars,
                    weighted = TRUE, tol = 1e-17, warn = TRUE, max.fail = 0,
                    verbose = FALSE, .ndone = 0)
          {
            if (missing(Nmif)) Nmif <- object@Nmif
            if (missing(start)) start <- coef(object)
            if (missing(pars)) pars <- object@pars
            if (missing(ivps)) ivps <- object@ivps
            if (missing(rw.sd)) rw.sd <- object@random.walk.sd
            if (missing(alg.pars)) alg.pars <- object@alg.pars
            theta <- start
            start.names <- names(start)
            if (is.null(start.names))
              stop("mif error: ",sQuote("start")," must be a named vector",call.=FALSE)

            sigma <- rep(0,length(start))
            names(sigma) <- start.names

            rw.names <- names(rw.sd)
            if (is.null(rw.names) || any(rw.sd<0))
              stop("mif error: ",sQuote("rw.sd")," must be a named non-negative numerical vector",call.=FALSE)
            if (!all(rw.names%in%start.names))
              stop("mif error: all the names of ",sQuote("rw.sd")," must be names of ",sQuote("start"),call.=FALSE)

            rw.names <- names(rw.sd[rw.sd>0])
            if (length(pars) == 0)
              stop("mif error: ",sQuote("pars")," must be a nonempty character vector",call.=FALSE)
            if (
                !is.character(pars) ||
                !is.character(ivps) ||
                !all(pars%in%start.names) ||
                !all(ivps%in%start.names) ||
                any(pars%in%ivps) ||
                any(ivps%in%pars) ||
                !all(pars%in%rw.names) ||
                !all(ivps%in%rw.names)
                )
              stop(
                   "mif error: ",
                   sQuote("pars")," and ",sQuote("ivps"),
                   " must be mutually disjoint subsets of ",
                   sQuote("names(start)"),
                   " and must have a positive random-walk SDs specified in ",
                   sQuote("rw.sd"),
                   call.=FALSE
                   )

            if (!all(rw.names%in%c(pars,ivps))) {
              extra.rws <- rw.names[!(rw.names%in%c(pars,ivps))]
              warning(
                      "mif warning: the variable(s) ",
                      paste(extra.rws,collapse=", "),
                      " have positive random-walk SDs specified, but are included in neither ",
                      sQuote("pars")," nor ",sQuote("ivps"),
                      ". These random walk SDs are ignored.",
                      call.=FALSE
                      )
            }
            rw.sd <- rw.sd[c(pars,ivps)]
            rw.names <- names(rw.sd)

            sigma[rw.names] <- rw.sd

            if (!all(c('Np','cooling.factor','ic.lag','var.factor')%in%names(alg.pars)))
              stop(
                   "mif error: ",sQuote("alg.pars")," must be a named list with elements ",
                   sQuote("Np"),",",sQuote("cooling.factor"),",",sQuote("ic.lag"),
                   ",and ",sQuote("var.factor"),
                   call.=FALSE
                   )
            
            conv.rec <- matrix(
                               data=NA,
                               nrow=Nmif+1,
                               ncol=length(theta)+2,
                               dimnames=list(
                                 seq(.ndone,.ndone+Nmif),
                                 c('loglik','nfail',names(theta))
                                 )
                               )
            conv.rec[1,] <- c(NA,NA,theta)
            
            for (n in seq(length=Nmif)) { # main loop

              ## compute the cooled sigma
              cool.sched <- try(
                                mif.cooling(alg.pars$cooling.factor,.ndone+n),
                                silent=FALSE
                                )
              if (inherits(cool.sched,'try-error'))
                stop("mif error: cooling schedule error",call.=FALSE)
              sigma.n <- sigma*cool.sched$alpha

              ## initialize the particles' parameter portion...
              P <- try(
                       particles(object,Np=alg.pars$Np,center=theta,sd=sigma.n*alg.pars$var.factor),
                       silent=FALSE
                       )
              if (inherits(P,'try-error'))
                stop("mif error: error in ",sQuote("particles"),call.=FALSE)

              ## run the particle filter
              x <- try(
                       pfilter.internal(
                                        object=object,
                                        params=P,
                                        tol=tol,
                                        warn=warn,
                                        max.fail=max.fail,
                                        pred.mean=(n==Nmif),
                                        pred.var=(weighted||(n==Nmif)),
                                        filter.mean=TRUE,
                                        .rw.sd=sigma.n[pars],
                                        verbose=verbose
                                        ),
                       silent=FALSE
                       )
              if (inherits(x,'try-error'))
                stop("mif error: error in ",sQuote("pfilter"),call.=FALSE)

              if (weighted) {           # MIF update rule
                v <- x$pred.var[pars,,drop=FALSE] # the prediction variance
                v1 <- cool.sched$gamma*(1+alg.pars$var.factor^2)*sigma[pars]^2
                theta.hat <- cbind(theta[pars],x$filter.mean[pars,,drop=FALSE])
                theta[pars] <- theta[pars]+apply(apply(theta.hat,1,diff)/t(v),2,sum)*v1
              } else {                  # unweighted (flat) average
                theta.hat <- x$filter.mean[pars,,drop=FALSE]
                theta[pars] <- apply(theta.hat,1,mean)
              }
              
              ## update the IVPs using fixed-lag smoothing
              theta[ivps] <- x$filter.mean[ivps,alg.pars$ic.lag]

              ## store a record of this iteration
              conv.rec[n+1,-c(1,2)] <- theta
              conv.rec[n,c(1,2)] <- c(x$loglik,x$nfail)

              if (verbose) cat("MIF iteration ",n," of ",Nmif," completed\n")

            }

            coef(object) <- theta

            new(
                "mif",
                object,
                ivps=ivps,
                pars=pars,
                Nmif=as.integer(Nmif),
                particles=object@particles,
                alg.pars=alg.pars,
                random.walk.sd=sigma[rw.names],
                pred.mean=x$pred.mean,
                pred.var=x$pred.var,
                filter.mean=x$filter.mean,
                conv.rec=conv.rec,
                cond.loglik=x$cond.loglik,
                eff.sample.size=x$eff.sample.size,
                loglik=x$loglik
                )
          }
          )

setMethod(
          'continue',
          'mif',
          function (object, Nmif, ...) {
            ndone <- object@Nmif
            obj <- mif(object,Nmif=Nmif,.ndone=ndone,...)
            object@conv.rec[ndone+1,c('loglik','nfail')] <- obj@conv.rec[1,c('loglik','nfail')]
            obj@conv.rec <- rbind(
                                  object@conv.rec,
                                  obj@conv.rec[-1,colnames(object@conv.rec)]
                                  )
            obj@Nmif <- as.integer(ndone+Nmif)
            obj
          }
          )

