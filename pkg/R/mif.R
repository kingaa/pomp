## MIF algorithm functions

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

mif.internal <- function (object, Nmif = 1,
                          start = NULL,
                          pars = NULL, ivps = NULL,
                          particles = NULL,
                          rw.sd = NULL,
                          Np = NULL, cooling.factor = NULL, var.factor = NULL, ic.lag = NULL, 
                          weighted = TRUE, tol = 1e-17, max.fail = 0,
                          verbose = FALSE, .ndone = 0) {
  is.mif <- is(object,"mif")
  if (is.null(start)) {
    start <- coef(object)
    if (length(start)==0)
      stop("mif error: ",sQuote("start")," must be specified if ",
           sQuote("coef(object)")," is NULL",
           call.=FALSE)
  } else if (length(start)==0)
    stop("mif error: ",sQuote("start")," must be a named vector",call.=FALSE)
  start.names <- names(start)
  if (is.null(start.names))
    stop("mif error: ",sQuote("start")," must be a named vector",call.=FALSE)

  if (is.null(rw.sd)) {
    if (is.mif)
      rw.sd <- object@random.walk.sd
    else
      stop("mif error: ",sQuote("rw.sd")," must be specified",call.=FALSE)
  }
  rw.names <- names(rw.sd)
  if (is.null(rw.names) || any(rw.sd<0))
    stop("mif error: ",sQuote("rw.sd")," must be a named non-negative numerical vector",call.=FALSE)
  if (!all(rw.names%in%start.names))
    stop("mif error: all the names of ",sQuote("rw.sd")," must be names of ",sQuote("start"),call.=FALSE)
  rw.names <- names(rw.sd[rw.sd>0])
  if (length(rw.names) == 0)
    stop("mif error: ",sQuote("rw.sd")," must have one positive entry for each parameter to be estimated",call.=FALSE)

  if (is.null(pars)) {
    if (is.mif) 
      pars <- object@pars
    else
      stop("mif error: ",sQuote("pars")," must be specified",call.=FALSE)
  }
  if (length(pars)==0)
    stop("mif error: at least one ordinary (non-IVP) parameter must be estimated",call.=FALSE)
  if (is.null(ivps)) {
    if (is.mif)
      ivps <- object@ivps
    else
      stop("mif error: ",sQuote("ivps")," must be specified",call.=FALSE)
  }      
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

  if (is.null(particles)) {
    if (is.mif) 
      particles <- object@particles
    else
      stop("mif error: ",sQuote("particles")," must be specified",call.=FALSE)
  }

  if (is.null(Np)) {
    if (is.mif) 
      Np <- object@alg.pars$Np
    else
      stop("mif error: ",sQuote("Np")," must be specified",call.=FALSE)
  }
  Np <- as.integer(Np)
  if ((length(Np)!=1)||(Np < 1))
    stop("mif error: ",sQuote("Np")," must be a positive integer",call.=FALSE)
  if (is.null(ic.lag)) {
    if (is.mif) 
      ic.lag <- object@alg.pars$ic.lag
    else
      stop("mif error: ",sQuote("ic.lag")," must be specified",call.=FALSE)
  }
  ic.lag <- as.integer(ic.lag)
  if ((length(ic.lag)!=1)||(ic.lag < 1))
    stop("mif error: ",sQuote("ic.lag")," must be a positive integer",call.=FALSE)
  if (is.null(cooling.factor)) {
    if (is.mif) 
      cooling.factor <- object@alg.pars$cooling.factor
    else
      stop("mif error: ",sQuote("cooling.factor")," must be specified",call.=FALSE)
  }
  if ((length(cooling.factor)!=1)||(cooling.factor < 0)||(cooling.factor>1))
    stop("mif error: ",sQuote("cooling.factor")," must be a number between 0 and 1",call.=FALSE)
  if (is.null(var.factor)) {
    if (is.mif) 
      var.factor <- object@alg.pars$var.factor
    else
      stop("mif error: ",sQuote("var.factor")," must be specified",call.=FALSE)
  }
  if ((length(var.factor)!=1)||(var.factor < 0))
    stop("mif error: ",sQuote("var.factor")," must be a positive number",call.=FALSE)

  Nmif <- as.integer(Nmif)
  if (Nmif<0)
    stop("mif error: ",sQuote("Nmif")," must be a positive integer",call.=FALSE)

  if (verbose) {
    cat("performing",Nmif,"MIF iteration(s) to estimate parameter(s)",
        paste(pars,collapse=", "))
    if (length(ivps)>0)
      cat(" and IVP(s)",paste(ivps,collapse=", "))
    cat(" using random-walk with SD\n")
    print(rw.sd)
    cat(
        "using",Np,"particles, variance factor",var.factor,
        "\ninitial condition smoothing lag",ic.lag,
        "and cooling factor",cooling.factor,"\n"
        )
  }

  theta <- start

  sigma <- rep(0,length(start))
  names(sigma) <- start.names

  rw.sd <- rw.sd[c(pars,ivps)]
  rw.names <- names(rw.sd)

  sigma[rw.names] <- rw.sd

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

  if (!all(is.finite(theta[c(pars,ivps)]))) {
    stop(
         sQuote("mif"),
         " error: cannot estimate non-finite parameters: ",
         paste(
               c(pars,ivps)[!is.finite(theta[c(pars,ivps)])],
               collapse=","
               ),
         call.=FALSE
         )
  }
  
  obj <- new(
             "mif",
             as(object,"pomp"),
             ivps=ivps,
             pars=pars,
             Nmif=0L,
             particles=particles,
             alg.pars=list(
               Np=Np,
               cooling.factor=cooling.factor,
               var.factor=var.factor,
               ic.lag=ic.lag
               ),
             random.walk.sd=sigma[rw.names],
             conv.rec=conv.rec
             )

  for (n in seq_len(Nmif)) { # main loop

    ## compute the cooled sigma
    cool.sched <- try(
                      mif.cooling(cooling.factor,.ndone+n),
                      silent=FALSE
                      )
    if (inherits(cool.sched,'try-error'))
      stop("mif error: cooling schedule error",call.=FALSE)
    sigma.n <- sigma*cool.sched$alpha

    ## initialize the particles' parameter portion...
    P <- try(
             particles(obj,Np=Np,center=theta,sd=sigma.n*var.factor),
             silent=FALSE
             )
    if (inherits(P,'try-error'))
      stop("mif error: error in ",sQuote("particles"),call.=FALSE)

    ## run the particle filter
    x <- try(
             pfilter.internal(
                              object=obj,
                              params=P,
                              tol=tol,
                              max.fail=max.fail,
                              pred.mean=(n==Nmif),
                              pred.var=(weighted||(n==Nmif)),
                              filter.mean=TRUE,
                              save.states=FALSE,
                              .rw.sd=sigma.n[pars],
                              verbose=verbose
                              ),
             silent=FALSE
             )
    if (inherits(x,'try-error'))
      stop("mif error: error in ",sQuote("pfilter"),call.=FALSE)

    if (weighted) {           # MIF update rule
      v <- x$pred.var[pars,,drop=FALSE] # the prediction variance
      v1 <- cool.sched$gamma*(1+var.factor^2)*sigma[pars]^2
      theta.hat <- cbind(theta[pars],x$filter.mean[pars,,drop=FALSE])
      theta[pars] <- theta[pars]+colSums(apply(theta.hat,1,diff)/t(v))*v1
    } else {                  # unweighted (flat) average
      theta.hat <- x$filter.mean[pars,,drop=FALSE]
      theta[pars] <- rowMeans(theta.hat)
    }
    
    ## update the IVPs using fixed-lag smoothing
    theta[ivps] <- x$filter.mean[ivps,ic.lag]

    ## store a record of this iteration
    conv.rec[n+1,-c(1,2)] <- theta
    conv.rec[n,c(1,2)] <- c(x$loglik,x$nfail)

    if (verbose) cat("MIF iteration ",n," of ",Nmif," completed\n")

  }

  coef(obj) <- theta

  if (Nmif>0) {
    obj@Nmif <- as.integer(Nmif)
    obj@conv.rec <- conv.rec
    obj@pred.mean <- x$pred.mean
    obj@pred.var <- x$pred.var
    obj@filter.mean <- x$filter.mean
    obj@eff.sample.size <- x$eff.sample.size
    obj@cond.loglik <- x$cond.loglik
    obj@loglik <- x$loglik
  }

  obj
}

mif <- function (object, ... )
  stop("function ",sQuote("mif")," is undefined for objects of class ",sQuote(class(object)))
setGeneric('mif')

continue <- function (object, ... )
  stop("function ",sQuote("continue")," is undefined for objects of class ",sQuote(class(object)))
setGeneric('continue')

setMethod(
          "mif",
          "pomp",
          function (object, Nmif = 1,
                    start,
                    pars, ivps = character(0),
                    particles, rw.sd,
                    Np, ic.lag, var.factor, cooling.factor,
                    weighted = TRUE, tol = 1e-17, max.fail = 0,
                    verbose = getOption("verbose"))
          {
            if (missing(start)) start <- NULL
            if (missing(rw.sd))
              stop("mif error: ",sQuote("rw.sd")," must be specified",call.=FALSE)
            if (missing(pars)) {
              rw.names <- names(rw.sd)[rw.sd>0]
              pars <- rw.names[!(rw.names%in%ivps)]
            }
            if (missing(Np))
              stop("mif error: ",sQuote("Np")," must be specified",call.=FALSE)
            if (missing(ic.lag))
              stop("mif error: ",sQuote("ic.lag")," must be specified",call.=FALSE)
            if (missing(var.factor))
              stop("mif error: ",sQuote("var.factor")," must be specified",call.=FALSE)
            if (missing(cooling.factor))
              stop("mif error: ",sQuote("cooling.factor")," must be specified",call.=FALSE)
              
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
              if (!all(c('Np','center','sd','...')%in%names(formals(particles))))
                stop(
                     "mif error: ",
                     sQuote("particles"),
                     " must be a function of prototype ",
                     sQuote("particles(Np,center,sd,...)"),
                     call.=FALSE
                     )
            }
              
            mif.internal(object,Nmif=Nmif,start=start,pars=pars,ivps=ivps,particles=particles,
                         rw.sd=rw.sd,Np=Np,cooling.factor=cooling.factor,
                         var.factor=var.factor,ic.lag=ic.lag,
                         weighted=weighted,tol=tol,max.fail=max.fail,
                         verbose=verbose,.ndone=0)

          }
          )
          

setMethod(
          "mif",
          "mif",
          function (object, Nmif, ...)
          {
            if (missing(Nmif)) Nmif <- object@Nmif
            mif.internal(object,Nmif=Nmif,...)
          }
          )

setMethod(
          'continue',
          'mif',
          function (object, Nmif = 1, ...) {
            ndone <- object@Nmif
            obj <- mif.internal(object,Nmif=Nmif,.ndone=ndone,...)
            object@conv.rec[ndone+1,c('loglik','nfail')] <- obj@conv.rec[1,c('loglik','nfail')]
            obj@conv.rec <- rbind(
                                  object@conv.rec,
                                  obj@conv.rec[-1,colnames(object@conv.rec)]
                                  )
            obj@Nmif <- as.integer(ndone+Nmif)
            obj
          }
          )
