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
                    pars = stop(sQuote("pars")," must be specified"),
                    ivps = character(0),
                    particles,
                    rw.sd = stop(sQuote("rw.sd")," must be specified"),
                    alg.pars = stop(sQuote("alg.pars")," must be specified"),
                    weighted = TRUE, tol = 1e-17, warn = TRUE, max.fail = 0,
                    verbose = FALSE, .ndone = 0) {
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
              stop("mif error: ",sQuote("particles")," must be a function of prototype ",sQuote("particles(Np,center,sd,...)"),call.=FALSE)
            if (missing(start)) {
              start <- coef(object)
              if (length(start)==0)
                stop("mif error: ",sQuote("start")," must be specified",call.=FALSE)
            }
            start.names <- names(start)
            if (is.null(start.names))
              stop("mif error: ",sQuote("start")," must be a named vector",call.=FALSE)
            if (length(pars) == 0)
              stop("mif error: ",sQuote("pars")," must be a nonempty character vector",call.=FALSE)
            if (
                !is.character(pars) ||
                !is.character(ivps) ||
                !all(pars%in%start.names) ||
                !all(ivps%in%start.names) ||
                any(pars%in%ivps) ||
                any(ivps%in%pars)
                )
              stop("mif error: ",sQuote("pars")," and ",sQuote("ivps")," must be mutually disjoint elements of ",sQuote("names(start)"),call.=FALSE)
            Nv <- length(start)
            if ((length(rw.sd)==1) && (rw.sd==0)) {
              rw.sd <- rep(0,Nv)
              names(rw.sd) <- start.names
            }
            rw.names <- names(rw.sd)
            if (any(!(rw.names%in%start.names)))
              stop("mif error: all the names of ",sQuote("rw.sd")," must be names of ",sQuote("start"),call.=FALSE)
            if (any(rw.sd[c(pars,ivps)]==0)) {
              zero.pars <- names(which(rw.sd[c(pars,ivps)]==0))
              stop(
                   "mif error: for every parameter you wish to estimate, you must specify a positive ",
                   sQuote("rw.sd"),
                   ".  ",
                   sQuote("rw.sd"),
                   " is zero for ",
                   paste(zero.pars,collapse=", "),
                   call.=FALSE
                   )
            }
            if (!all(c('Np','cooling.factor','ic.lag','var.factor')%in%names(alg.pars)))
              stop(
                   "mif error: ",sQuote("alg.pars")," must be a named list with elements ",sQuote("Np"),",",sQuote("cooling.factor"),",",sQuote("ic.lag"),",and ",sQuote("var.factor"),
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
                              c('loglik','nfail',names(start))
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
          function (object, Nmif = object@Nmif, start = coef(object),
                    pars = object@pars, ivps = object@ivps, rw.sd = object@random.walk.sd,
                    alg.pars = object@alg.pars,
                    weighted = TRUE, tol = 1e-17, warn = TRUE, max.fail = 0,
                    verbose = FALSE, .ndone = 0) {
            theta <- start
            sigma <- rep(0,length(start))
            names(sigma) <- names(start)
            rw.names <- names(rw.sd)
            if (!all(rw.names%in%names(start)))
              stop("mif error: all the names of ",sQuote("rw.sd")," must be names of ",sQuote("start"),call.=FALSE)
            sigma[rw.names] <- rw.sd
            if (!all(c('Np','cooling.factor','ic.lag','var.factor')%in%names(alg.pars)))
              stop("mif error: ",sQuote("alg.pars")," must be a named list with elements ",sQuote("Np"),",",sQuote("cooling.factor"),",",sQuote("ic.lag"),",and ",sQuote("var.factor"),call.=FALSE)
            conv.rec <- matrix(NA,
                               nrow=Nmif+1,
                               ncol=length(theta)+2,
                               dimnames=list(
                                 seq(.ndone,.ndone+Nmif),
                                 c('loglik','nfail',names(theta))
                                 )
                               )
            conv.rec[1,] <- c(NA,NA,theta)
            for (n in seq(length=Nmif)) {
              cool.sched <- try(
                                mif.cooling(alg.pars$cooling.factor,.ndone+n),
                                silent=FALSE
                                )
              if (inherits(cool.sched,'try-error'))
                stop("mif error: cooling schedule error",call.=FALSE)
              sigma.n <- sigma*cool.sched$alpha
              P <- try(
                       particles(object,Np=alg.pars$Np,center=theta,sd=sigma.n*alg.pars$var.factor),
                       silent=FALSE
                       )
              if (inherits(P,'try-error'))
                stop("mif error: error in ",sQuote("particles"),call.=FALSE)
              X <- try(
                       init.state(object,params=P),
                       silent=FALSE
                       )
              if (inherits(X,'try-error'))
                stop("mif error: error in ",sQuote("init.state"),call.=FALSE)
              x <- try(
                       pfilter(
                               as(object,'pomp'),
                               xstart=X,
                               params=P,
                               tol=tol,
                               warn=warn,
                               max.fail=max.fail,
                               pred.mean=(n==Nmif),
                               pred.var=TRUE,
                               filter.mean=TRUE,
                               .rw.sd=sigma.n[pars]
                               ),
                       silent=FALSE
                       )
              if (inherits(x,'try-error'))
                stop("mif error: error in ",sQuote("pfilter"),call.=FALSE)

              v <- x$pred.var[pars,,drop=FALSE]
              
              if (weighted) {                     # MIF update rule
                v1 <- cool.sched$gamma*(1+alg.pars$var.factor^2)*sigma[pars]^2
                theta.hat <- cbind(theta[pars],x$filter.mean[pars,,drop=FALSE])
                theta[pars] <- theta[pars]+apply(apply(theta.hat,1,diff)/t(v),2,sum)*v1
              } else {                            # unweighted (flat) average
                theta.hat <- x$filter.mean[pars,,drop=FALSE]
                theta[pars] <- apply(theta.hat,1,mean)
              }
              theta[ivps] <- x$filter.mean[ivps,alg.pars$ic.lag]
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
                random.walk.sd=sigma,
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
            object@conv.rec[ndone+1,'loglik'] <- obj@conv.rec[1,'loglik']
            obj@conv.rec <- rbind(
                                  object@conv.rec,
                                  obj@conv.rec[-1,colnames(object@conv.rec)]
                                  )
            obj@Nmif <- as.integer(ndone+Nmif)
            obj
          }
          )

