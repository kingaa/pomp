## MIF algorithm functions

default.pomp.particles.fun <- function (Np, center, sd, ...) {
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

cooling.function <- function (type, perobs, fraction, ntimes) {
  switch(
         type,
         geometric={
           factor <- fraction^(1/50)
           if (perobs) {
             function (nt, m) {
               alpha <- factor^(nt/ntimes+m-1)
               list(alpha=alpha,gamma=alpha^2)
             }
           } else {
             function (nt, m) {
               alpha <- factor^(m-1)
               list(alpha=alpha,gamma=alpha^2)
             }
           }
         },
         hyperbolic={
           if (perobs) {
             scal <- (50*ntimes*fraction-1)/(1-fraction)
             function (nt, m) {
               alpha <- (1+scal)/(scal+nt+ntimes*(m-1))
               list(alpha=alpha,gamma=alpha^2)
             }
           } else {
             scal <- (50*fraction-1)/(1-fraction)
             function (nt, m) {
               alpha <- (1+scal)/(scal+m-1)
               list(alpha=alpha,gamma=alpha^2)
             }

           }
         },
         stop("unrecognized cooling schedule type ",sQuote(type))
         )
}

mif.cooling <- function (factor, n) {   # default geometric cooling schedule
  alpha <- factor^(n-1)
  list(alpha=alpha,gamma=alpha^2)
}

mif2.cooling <- function (frac, nt, m, n) {   # cooling schedule for mif2
  ## frac is the fraction of cooling after 50 iterations
  scal <- (50*n*frac-1)/(1-frac)
  alpha <- (1+scal)/(scal+nt+n*(m-1))
  list(alpha=alpha)
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

mif.internal <- function (object, Nmif,
                          start, pars, ivps,
                          particles,
                          rw.sd,
                          Np, var.factor, ic.lag,
                          cooling.type, cooling.fraction, cooling.factor, 
                          method,
                          tol, max.fail,
                          verbose, transform, .ndone = 0,
                          paramMatrix) {
  
  transform <- as.logical(transform)
  
  if (length(start)==0)
    stop(
         "mif error: ",sQuote("start")," must be specified if ",
         sQuote("coef(object)")," is NULL",
         call.=FALSE
         )
  
  if (transform)
    start <- partrans(object,start,dir="inverse")
  
  start.names <- names(start)
  if (is.null(start.names))
    stop("mif error: ",sQuote("start")," must be a named vector",call.=FALSE)
  
  if (missing(rw.sd))
    stop("mif error: ",sQuote("rw.sd")," must be specified",call.=FALSE)
  
  rw.names <- names(rw.sd)
  if (is.null(rw.names) || any(rw.sd<0))
    stop("mif error: ",sQuote("rw.sd")," must be a named non-negative numerical vector",call.=FALSE)
  if (!all(rw.names%in%start.names))
    stop("mif error: all the names of ",sQuote("rw.sd")," must be names of ",sQuote("start"),call.=FALSE)
  rw.names <- names(rw.sd[rw.sd>0])
  if (length(rw.names) == 0)
    stop("mif error: ",sQuote("rw.sd")," must have one positive entry for each parameter to be estimated",call.=FALSE)
  
  if (missing(pars)) stop("mif error: ",sQuote("pars")," must be specified",call.=FALSE)
  if (missing(ivps)) stop("mif error: ",sQuote("ivps")," must be specified",call.=FALSE)
  
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
            ngettext(length(extra.rws),"mif warning: the variable ",
                     "mif warning: the variables "),
            paste(sQuote(extra.rws),collapse=", "),
            ngettext(length(extra.rws)," has positive random-walk SD specified, but is included in neither ",
                     " have positive random-walk SDs specified, but are included in neither "),
            sQuote("pars")," nor ",sQuote("ivps"),
            ngettext(length(extra.rws),". This random walk SD will be ignored.",
                     ". These random walk SDs will be ignored."),
            call.=FALSE
            )
  }
  rw.sd <- rw.sd[c(pars,ivps)]
  rw.names <- names(rw.sd)
  
  if (missing(particles))
    stop("mif error: ",sQuote("particles")," must be specified",call.=FALSE)
  
  ntimes <- length(time(object))
  if (missing(Np)) stop("mif error: ",sQuote("Np")," must be specified",call.=FALSE)
  if (is.function(Np)) {
    Np <- try(
              vapply(seq.int(from=0,to=ntimes,by=1),Np,numeric(1)),
              silent=FALSE
              )
    if (inherits(Np,"try-error"))
      stop("if ",sQuote("Np")," is a function, it must return a single positive integer")
  }
  if (length(Np)==1)
    Np <- rep(Np,times=ntimes+1)
  else if (length(Np)!=(ntimes+1))
    stop(sQuote("Np")," must have length 1 or length ",ntimes+1)
  if (any(Np<=0))
    stop("number of particles, ",sQuote("Np"),", must always be positive")
  if (!is.numeric(Np))
    stop(sQuote("Np")," must be a number, a vector of numbers, or a function")
  Np <- as.integer(Np)
  
  if (missing(ic.lag)) stop("mif error: ",sQuote("ic.lag")," must be specified",call.=FALSE)
  ic.lag <- as.integer(ic.lag)
  if ((length(ic.lag)!=1)||(ic.lag<1))
    stop("mif error: ",sQuote("ic.lag")," must be a positive integer",call.=FALSE)
  if (ic.lag>ntimes) {
    warning(
            "mif warning: ",sQuote("ic.lag")," = ",ic.lag," > ",ntimes,
            " = length(time(",sQuote("object"),"))",
            " is nonsensical.  Setting ",sQuote("ic.lag")," = ",ntimes,".",
            call.=FALSE
            )
    ic.lag <- length(time(object))
  }
  if ((length(pars)==0)&&(ic.lag<length(time(object)))) {
    warning(
            "mif warning: only IVPs are to be estimated, yet ",sQuote("ic.lag")," = ",ic.lag,
            " < ",ntimes," = length(time(",sQuote("object"),")),",
            " so unnecessary work is to be done.",
            call.=FALSE
            )
  }
  
  ## the following deals with the deprecated option 'cooling.factor'
  if (!missing(cooling.factor)) {
    warning(sQuote("cooling.factor")," is deprecated.\n",
            "See ?mif for instructions on specifying the cooling schedule.",
            call.=FALSE)
    cooling.factor <- as.numeric(cooling.factor)
    if ((length(cooling.factor)!=1)||(cooling.factor<0)||(cooling.factor>1))
      stop("mif error: ",sQuote("cooling.factor")," must be a number between 0 and 1",call.=FALSE)
    if (missing(cooling.fraction)) {
      cooling.fraction <- cooling.factor^50
    } else {
      warning("specification of ",sQuote("cooling.factor"),
              " is overridden by that of ",sQuote("cooling.fraction"),
              call.=FALSE)
    }
  }

  if (missing(cooling.fraction))
    stop("mif error: ",sQuote("cooling.fraction")," must be specified",call.=FALSE)
  cooling.fraction <- as.numeric(cooling.fraction)
  if ((length(cooling.fraction)!=1)||(cooling.fraction<0)||(cooling.fraction>1))
    stop("mif error: ",sQuote("cooling.fraction")," must be a number between 0 and 1",call.=FALSE)
  
  cooling <- cooling.function(
                              type=cooling.type,
                              perobs=(method=="mif2"),
                              fraction=cooling.fraction,
                              ntimes=ntimes
                              )

  if ((method=="mif2")&&(Np[1L]!=Np[ntimes+1]))
    stop("the first and last values of ",sQuote("Np")," must agree when method = ",sQuote("mif2"))
  
  if (missing(var.factor))
    stop("mif error: ",sQuote("var.factor")," must be specified",call.=FALSE)
  if ((length(var.factor)!=1)||(var.factor < 0))
    stop("mif error: ",sQuote("var.factor")," must be a positive number",call.=FALSE)
  
  if (missing(Nmif))
    stop("mif error: ",sQuote("Nmif")," must be specified",call.=FALSE)
  Nmif <- as.integer(Nmif)
  if (Nmif<0)
    stop("mif error: ",sQuote("Nmif")," must be a positive integer",call.=FALSE)

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
  conv.rec[1L,] <- c(NA,NA,theta)
  
  if (!all(is.finite(theta[c(pars,ivps)]))) {
    stop(
         sQuote("mif")," cannot estimate non-finite parameters.\n",
         "The following ",if (transform) "transformed ", "parameters are non-finite: ",
         paste(
               sQuote(c(pars,ivps)[!is.finite(theta[c(pars,ivps)])]),
               collapse=","
               ),
         call.=FALSE
         )
  }
  
  obj <- as(object,"pomp")
  
  if (Nmif>0) {
    tmp.mif <- new("mif",object,particles=particles,Np=Np[1L])
  } else {
    pfp <- obj
  }
  
  have.parmat <- !(missing(paramMatrix) || length(paramMatrix)==0)

  for (n in seq_len(Nmif)) { ## iterate the filtering

    ## get the intensity of artificial noise from the cooling schedule
    cool.sched <- cooling(nt=1,m=.ndone+n)
    sigma.n <- sigma*cool.sched$alpha
    
    ## initialize the parameter portions of the particles
    P <- try(
             particles(
                       tmp.mif,
                       Np=Np[1L],
                       center=theta,
                       sd=sigma.n*var.factor
                       ),
             silent = FALSE
             )
    if (inherits(P,"try-error")) 
      stop("mif error: error in ",sQuote("particles"),call.=FALSE)

    if ((method=="mif2") && ((n>1) || have.parmat)) {
      ## use pre-existing particle matrix
      P[pars,] <- paramMatrix[pars,]
    }

    pfp <- try(
               pfilter.internal(
                                object=obj,
                                params=P, 
                                Np=Np,
                                tol=tol,
                                max.fail=max.fail,
                                pred.mean=(n==Nmif),
                                pred.var=((method=="mif")||(n==Nmif)),
                                filter.mean=TRUE,
                                cooling=cooling,
                                cooling.m=.ndone+n,
                                .mif2=(method=="mif2"),
                                .rw.sd=sigma.n[pars],
                                .transform=transform,
                                save.states=FALSE, 
                                save.params=FALSE,
                                verbose=verbose
                                ),
               silent=FALSE
               )
    if (inherits(pfp,"try-error")) 
      stop("mif error: error in ",sQuote("pfilter"),call.=FALSE)

    ## update parameters
    switch(
           method,
           mif={              # original Ionides et al. (2006) average
             theta <- .Call(mif_update,pfp,theta,cool.sched$gamma,var.factor,sigma,pars)
           },
           unweighted={                 # unweighted average
             theta[pars] <- rowMeans(pfp@filter.mean[pars,,drop=FALSE])
           },
           fp={                         # fixed-point iteration
             theta[pars] <- pfp@filter.mean[pars,ntimes,drop=FALSE]
           },
           mif2={                     # "efficient" iterated filtering
             paramMatrix <- pfp@paramMatrix
             theta[pars] <- rowMeans(paramMatrix[pars,,drop=FALSE])
           },
           stop("unrecognized method ",sQuote(method))
           )
    theta[ivps] <- pfp@filter.mean[ivps,ic.lag]
    conv.rec[n+1,-c(1,2)] <- theta
    conv.rec[n,c(1,2)] <- c(pfp@loglik,pfp@nfail)
    
    if (verbose) cat("MIF iteration ",n," of ",Nmif," completed\n")
    
  } ### end of main loop

  ## back transform the parameter estimate if necessary
  if (transform) theta <- partrans(pfp,theta,dir="forward")
  
  new(
      "mif",
      pfp,
      transform=transform,
      params=theta,
      ivps=ivps,
      pars=pars,
      Nmif=Nmif,
      particles=particles,
      var.factor=var.factor,
      ic.lag=ic.lag,
      random.walk.sd=sigma[rw.names],
      tol=tol,
      conv.rec=conv.rec,
      method=method,
      cooling.type=cooling.type,
      cooling.fraction=cooling.fraction,
      paramMatrix=if (method=="mif2") paramMatrix else array(data=numeric(0),dim=c(0,0))
      )
}

setGeneric('mif',function(object,...)standardGeneric("mif"))

setMethod(
          "mif",
          signature=signature(object="pomp"),
          function (object, Nmif = 1,
                    start,
                    pars, ivps = character(0),
                    particles, rw.sd,
                    Np, ic.lag, var.factor,
                    cooling.type = c("geometric","hyperbolic"),
                    cooling.fraction, cooling.factor,
                    method = c("mif","unweighted","fp","mif2"),
                    tol = 1e-17, max.fail = Inf,
                    verbose = getOption("verbose"),
                    transform = FALSE,
                    ...) {
            
            transform <- as.logical(transform)
            method <- match.arg(method)
            
            if (missing(start)) start <- coef(object)
            if (missing(rw.sd))
              stop("mif error: ",sQuote("rw.sd")," must be specified",call.=FALSE)
            if (missing(pars)) {
              rw.names <- names(rw.sd)[rw.sd>0]
              pars <- rw.names[!(rw.names%in%ivps)]
            }
            if (missing(Np))
              stop("mif error: ",sQuote("Np")," must be specified",call.=FALSE)
            if (missing(ic.lag) && length(ivps)>0)
              stop("mif error: ",sQuote("ic.lag")," must be specified if ",sQuote("ivps")," are",call.=FALSE)
            if (missing(var.factor))
              stop("mif error: ",sQuote("var.factor")," must be specified",call.=FALSE)

            cooling.type <- match.arg(cooling.type)
            
            if (missing(particles)) { # use default: normal distribution
              particles <- default.pomp.particles.fun
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
            
            mif.internal(
                         object=object,
                         Nmif=Nmif,
                         start=start,
                         pars=pars,
                         ivps=ivps,
                         particles=particles,
                         rw.sd=rw.sd,
                         Np=Np,
                         cooling.type=cooling.type,
                         cooling.factor=cooling.factor,
                         cooling.fraction=cooling.fraction,
                         var.factor=var.factor,
                         ic.lag=ic.lag,
                         method=method,
                         tol=tol,
                         max.fail=max.fail,
                         verbose=verbose,
                         transform=transform
                         )
            
          }
          )


setMethod(
          "mif",
          signature=signature(object="pfilterd.pomp"),
          function (object, Nmif = 1, Np, tol,
                    ...) {
            
            if (missing(Np)) Np <- object@Np
            if (missing(tol)) tol <- object@tol
            
            mif(
                object=as(object,"pomp"),
                Nmif=Nmif,
                Np=Np,
                tol=tol,
                ...
                )
          }
          )

setMethod(
          "mif",
          signature=signature(object="mif"),
          function (object, Nmif,
                    start,
                    pars, ivps,
                    particles, rw.sd,
                    Np, ic.lag, var.factor,
                    cooling.type, cooling.fraction,
                    method,
                    tol,
                    transform,
                    ...) {
            
            if (missing(Nmif)) Nmif <- object@Nmif
            if (missing(start)) start <- coef(object)
            if (missing(pars)) pars <- object@pars
            if (missing(ivps)) ivps <- object@ivps
            if (missing(particles)) particles <- object@particles
            if (missing(rw.sd)) rw.sd <- object@random.walk.sd
            if (missing(ic.lag)) ic.lag <- object@ic.lag
            if (missing(var.factor)) var.factor <- object@var.factor
            if (missing(cooling.type)) cooling.type <- object@cooling.type
            if (missing(cooling.fraction)) cooling.fraction <- object@cooling.fraction
            if (missing(method)) method <- object@method
            if (missing(transform)) transform <- object@transform
            transform <- as.logical(transform)

            if (missing(Np)) Np <- object@Np
            if (missing(tol)) tol <- object@tol

            mif(
                object=as(object,"pomp"),
                Nmif=Nmif,
                start=start,
                pars=pars,
                ivps=ivps,
                particles=particles,
                rw.sd=rw.sd,
                Np=Np,
                cooling.type=cooling.type,
                cooling.fraction=cooling.fraction,
                var.factor=var.factor,
                ic.lag=ic.lag,
                method=method,
                tol=tol,
                transform=transform,
                ...
                )
          }
          )

setMethod(
          'continue',
          signature=signature(object='mif'),
          function (object, Nmif = 1,
                    ...) {
            
            ndone <- object@Nmif
            
            obj <- mif(
                       object=object,
                       Nmif=Nmif,
                       .ndone=ndone,
                       paramMatrix=object@paramMatrix,
                       ...
                       )
            
            object@conv.rec[ndone+1,c('loglik','nfail')] <- obj@conv.rec[1L,c('loglik','nfail')]
            obj@conv.rec <- rbind(
                                  object@conv.rec,
                                  obj@conv.rec[-1L,colnames(object@conv.rec)]
                                  )
            obj@Nmif <- as.integer(ndone+Nmif)
            
            obj
          }
          )

mif.profileDesign <- function (object, profile, lower, upper, nprof, ivps, 
                               rw.sd, Np, ic.lag, var.factor, cooling.factor,option, cooling.fraction, paramMatrix, ...)
{
  if (missing(profile)) profile <- list()
  if (missing(lower)) lower <- numeric(0)
  if (missing(upper)) upper <- lower
  if (length(lower)!=length(upper))
    stop(sQuote("lower")," and ",sQuote("upper")," must be of the same length")
  pars <- names(lower)
  if (missing(ivps)) ivps <- character(0)
  Np <- as.integer(Np)
  
  pd <- do.call(profileDesign,c(profile,list(lower=lower,upper=upper,nprof=nprof)))
  
  object <- as(object,"pomp")
  
  pp <- coef(object)
  idx <- !(names(pp)%in%names(pd))
  if (any(idx)) pd <- cbind(pd,as.list(pp[idx]))
  
  ans <- vector(mode="list",length=nrow(pd))
  for (k in seq_len(nrow(pd))) {
    ans[[k]] <- list(
                     mf=mif(
                       object,
                       Nmif=0,
                       start=unlist(pd[k,]),
                       pars=pars,
                       ivps=ivps,
                       rw.sd=rw.sd,
                       Np=Np,
                       ic.lag=ic.lag,
                       var.factor=var.factor,
                       cooling.factor=cooling.factor,
                       option=option,
                       cooling.fraction=cooling.fraction,
                       paramMatrix=paramMatrix,
                       ...
                       )
                     )
  }
  
  ans
}
