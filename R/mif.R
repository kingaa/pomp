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

mif.pomp <- function (object, Nmif = 1,
                      start,
                      pars = stop("'pars' must be specified"),
                      ivps = character(0),
                      particles,
                      rw.sd = stop("'rw.sd' must be specified"),
                      alg.pars = stop("'alg.pars' must be specified"),
                      weighted = TRUE, tol = 1e-17, warn = TRUE, max.fail = 0,
                      .ndone = 0) {
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
    stop("'particles' must be a function of prototype 'particles(Np,center,sd,...)'")
  if (missing(start)) {
    if (length(coef(object))>0) {
      start <- coef(object)
    } else {
      stop("'start' must be specified")
    }
  }
  start.names <- names(start)
  if (is.null(start.names))
    stop("mif error: 'start' must be a named vector")
  if (length(pars) == 0)
    stop("mif error: 'pars' must be a nonempty character vector")
  if (
      !is.character(pars) ||
      !is.character(ivps) ||
      !all(pars%in%start.names) ||
      !all(ivps%in%start.names) ||
      any(pars%in%ivps) ||
      any(ivps%in%pars)
      )
    stop("'pars' and 'ivps' must be mutually disjoint elements of 'names(start)'")
  Nv <- length(start)
  if ((length(rw.sd)==1) && (rw.sd==0)) {
    rw.sd <- rep(0,Nv)
    names(rw.sd) <- start.names
  }
  rw.names <- names(rw.sd)
  if (any(!(rw.names%in%start.names)))
    stop("all the names of 'rw.sd' must be names of 'start'")
  if (!all(c('Np','cooling.factor','ic.lag','var.factor')%in%names(alg.pars)))
    stop("'alg.pars' must be a named list with elements 'Np','cooling.factor','ic.lag',and 'var.factor'")
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
    mif.mif(newmif,Nmif=Nmif,weighted=weighted,tol=tol,warn=warn,max.fail=max.fail,.ndone=0)
  } else {
    newmif
  }
}
  
mif.mif <- function (object, Nmif = object@Nmif, start = coef(object),
                     pars = object@pars, ivps = object@ivps, rw.sd = object@random.walk.sd,
                     alg.pars = object@alg.pars,
                     weighted = TRUE, tol = 1e-17, warn = TRUE, max.fail = 0,
                     .ndone = 0) {
  theta <- start
  sigma <- rep(0,length(start))
  names(sigma) <- names(start)
  rw.names <- names(rw.sd)
  if (!all(rw.names%in%names(start)))
    stop("all the names of 'rw.sd' must be names of 'start'")
  sigma[rw.names] <- rw.sd
  if (!all(c('Np','cooling.factor','ic.lag','var.factor')%in%names(alg.pars)))
    stop("'alg.pars' must be a named list with elements 'Np','cooling.factor','ic.lag',and 'var.factor'")
  conv.rec <- matrix(NA,
                     nrow=Nmif+1,
                     ncol=length(theta)+2,
                     dimnames=list(
                       seq(0,Nmif),
                       c('loglik','nfail',names(theta))
                       )
                     )
  conv.rec[1,] <- c(NA,NA,theta)
  for (n in seq(length=Nmif)) {
    cool.sched <- try(
                      mif.cooling(alg.pars$cooling.factor,.ndone+n),
                      silent=T
                      )
    if (inherits(cool.sched,'try-error'))
      stop("mif error: cooling schedule error\n",cool.sched)
    sigma.n <- sigma*cool.sched$alpha
    P <- try(
             particles(object,Np=alg.pars$Np,center=theta,sd=sigma.n*alg.pars$var.factor),
             silent=FALSE
             )
    if (inherits(P,'try-error'))
      stop("mif error: error in 'particles'")
    X <- try(
             init.state(object,params=P,t0=object@t0),
             silent=FALSE
             )
    if (inherits(X,'try-error'))
      stop("mif error: error in 'init.state'")
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
      stop("mif error: error in 'pfilter'")

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
    conv.rec[n+1,-1] <- c(x$nfail,theta)
    conv.rec[n,1] <- x$loglik
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

setMethod('mif','pomp',mif.pomp)
setMethod('mif','mif',mif.mif)

setMethod(
          'continue',
          'mif',
          function (object, Nmif, ...) {
            ndone <- object@Nmif
            obj <- mif.mif(object,Nmif=Nmif,.ndone=ndone,...)
            object@conv.rec[ndone+1,'loglik'] <- obj@conv.rec[1,'loglik']
            obj@conv.rec <- rbind(
                                  object@conv.rec,
                                  obj@conv.rec[-1,colnames(object@conv.rec)]
                                  )
            obj@Nmif <- as.integer(ndone+Nmif)
            obj
          }
          )

