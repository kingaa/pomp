## MIF2 algorithm functions

## define a class of perturbation magnitudes
setClass(
         "mif2.perturb.sd",
         slots=c(
           sds="list",
           rwnames="character"
           ),
         prototype=prototype(
           sds=list(),
           rwnames=character(0)
           )
         )

## define the mif2d.pomp class
setClass(
         'mif2d.pomp',
         contains='pfilterd.pomp',
         slots=c(
           Nmif = 'integer',
           transform = 'logical',
           perturb.fn = 'function',
           rw.sd = 'mif2.perturb.sd',
           conv.rec = 'matrix'
           )
         )

mif2.sd <- function (...) {
  sds <- list(...)
  if (length(sds)==0)
    stop(sQuote("mif2.sd")," error: at least one parameter should be perturbed",call.=FALSE)
  rwnames <- names(sds)
  if (is.null(rwnames) || any(names(rwnames)==""))
    stop(sQuote("mif2.sd")," accepts only named arguments",call.=FALSE)
  for (n in rwnames) {
    sds[[n]] <- try(match.fun(sds[[n]]),silent=TRUE)
    if (inherits(sds[[n]],"try-error"))
      stop(sQuote("mif2.sd")," error: ",sQuote(n)," is not a function",call.=FALSE)
  }
  new("mif2.perturb.sd",sds=sds,rwnames=rwnames)
}

setGeneric("cooling_fn",function(object,...)standardGeneric("cooling_fn"))

setMethod("cooling_fn",
          signature=signature(object="mif2.perturb.sd"),
          definition=function (object, paramnames) {
            if (!all(object@rwnames %in% paramnames)) {
              unrec <- object@rwnames[!object@rwnames %in% paramnames]
              stop(sQuote("mif2")," error: the following parameter(s), ",
                   "which are supposed to be estimated, are not present: ",
                   paste(sapply(sQuote,unrec),collapse=","),
                   call.=FALSE)
            }
            function (mifiter, timept) {
              rw.sd <- setNames(numeric(length(object@rwnames)),object@rwnames)
              for (nm in object@rwnames) {
                rw.sd[nm] <- object@sds[[nm]](mifiter,timept)
              }
              rw.sd
            }
          })

geomcool <- function (sd, cooling.fraction.50 = 1) {
  if (missing(sd))
    stop(sQuote("geomcool")," error: ",sQuote("sd")," must be supplied",call.=FALSE)
  sd <- as.numeric(sd)
  if (sd <= 0)
    stop(sQuote("geomcool")," error: ",sQuote("sd")," must be non-negative",call.=FALSE)
  cooling.fraction.50 <- as.numeric(cooling.fraction.50)
  if (cooling.fraction.50 <= 0 || cooling.fraction.50 > 1)
    stop(sQuote("cooling.fraction.50")," must be in (0,1]",call.=FALSE)
  factor <- cooling.fraction.50^(1/50)
  function (mifiter, timept) {
    sd*(factor^(mifiter-1))
  }
}

hypcool <- function (sd, cooling.fraction.50 = 1, ntimes = NULL) {
  if (missing(sd))
    stop(sQuote("hypcool")," error: ",sQuote("sd")," must be supplied",call.=FALSE)
  sd <- as.numeric(sd)
  if (sd <= 0)
    stop(sQuote("hypcool")," error: ",sQuote("sd")," must be non-negative",call.=FALSE)
  cooling.fraction.50 <- as.numeric(cooling.fraction.50)
  if (cooling.fraction.50 <= 0 || cooling.fraction.50 > 1)
    stop(sQuote("cooling.fraction.50")," must be in (0,1]",call.=FALSE)
  if (is.null(ntimes)) {
    scal <- (50*cooling.fraction.50-1)/(1-cooling.fraction.50)
    function (mifiter, timept) {
      (1+scal)/(mifiter+scal)
    }
  } else {
    scal <- (50*ntimes*cooling.fraction.50-1)/(1-cooling.fraction.50)
    function (mifiter, timept) {
      sd*(1+scal)/(timept+ntimes*(mifiter-1)+scal)
    }
  }
}

ivphypcool <- function (sd, cooling.fraction.50 = 1) {
  if (missing(sd))
    stop(sQuote("ivphypcool")," error: ",sQuote("sd")," must be supplied",call.=FALSE)
  sd <- as.numeric(sd)
  if (sd <= 0)
    stop(sQuote("ivphypcool")," error: ",sQuote("sd")," must be non-negative",call.=FALSE)
  cooling.fraction.50 <- as.numeric(cooling.fraction.50)
  if (cooling.fraction.50 <= 0 || cooling.fraction.50 > 1)
    stop(sQuote("cooling.fraction.50")," must be in (0,1]",call.=FALSE)
  scal <- (50*cooling.fraction.50-1)/(1-cooling.fraction.50)
  function (mifiter, timept) {
    if (timept==1L)
      sd*(1+scal)/(mifiter+scal)
    else 0.0
  }
}

ivpgeomcool <- function (sd, cooling.fraction.50 = 1) {
  if (missing(sd))
    stop(sQuote("ivpgeomcool")," error: ",sQuote("sd")," must be supplied",call.=FALSE)
  sd <- as.numeric(sd)
  if (sd <= 0)
    stop(sQuote("ivpgeomcool")," error: ",sQuote("sd")," must be non-negative",call.=FALSE)
  cooling.fraction.50 <- as.numeric(cooling.fraction.50)
  if (cooling.fraction.50 <= 0 || cooling.fraction.50 > 1)
    stop(sQuote("cooling.fraction.50")," must be in (0,1]",call.=FALSE)
  factor <- cooling.fraction.50^(1/50)
  function (mifiter, timept) {
    if (timept==1L)
      sd*factor^(mifiter-1)
    else 0.0
  }
}

mif2.pfilter <- function (object, params, Np,
                          mifiter, cooling.fn, perturb.fn,
                          tol = 1e-17, max.fail = Inf,
                          transform = FALSE, verbose = FALSE,
                          filter.mean = FALSE,
                          .getnativesymbolinfo = TRUE) {

  gnsi <- as.logical(.getnativesymbolinfo)
  transform <- as.logical(transform)
  verbose <- as.logical(verbose)
  filter.mean <- as.logical(filter.mean)
  mifiter <- as.integer(mifiter)
  Np <- as.integer(Np)

  times <- time(object,t0=TRUE)
  ntimes <- length(times)-1

  paramnames <- rownames(params)
  npars <- nrow(params)

  loglik <- rep(NA,ntimes)
  eff.sample.size <- numeric(ntimes)
  nfail <- 0

  for (nt in seq_len(ntimes)) {

    ## perturb parameters
    params <- perturb.fn(params,sd=cooling.fn(mifiter,nt))

    if (transform)
      tparams <- partrans(object,params,dir="fromEstimationScale",
                          .getnativesymbolinfo=gnsi)

    if (nt == 1L) {
      ## get initial states
      x <- init.state(object,params=if (transform) tparams else params)

      if (filter.mean)
        filt.m <- array(dim=c(nrow(x),ntimes),
                        dimnames=list(rownames(x),NULL))
      else
        filt.m <- array(dim=c(0,0))
    }

    ## advance the state variables according to the process model
    X <- try(
             rprocess(
                      object,
                      xstart=x,
                      times=times[c(nt,nt+1)],
                      params=if (transform) tparams else params,
                      offset=1,
                      .getnativesymbolinfo=gnsi
                      ),
             silent=FALSE
             )
    if (inherits(X,'try-error'))
      stop(sQuote("mif2.pfilter")," error: process simulation error")

    ## determine the weights
    weights <- try(
                   dmeasure(
                            object,
                            y=object@data[,nt,drop=FALSE],
                            x=X,
                            times=times[nt+1],
                            params=if (transform) tparams else params,
                            log=FALSE,
                            .getnativesymbolinfo=gnsi
                            ),
                   silent=FALSE
                   )
    if (inherits(weights,'try-error'))
      stop(sQuote("mif2.pfilter")," error: error in calculation of weights",call.=FALSE)
    if (any(!is.finite(weights))) {
      stop(sQuote("mif2.pfilter")," error: ",sQuote("dmeasure"),
           " returns non-finite value",call.=FALSE)
    }
    gnsi <- FALSE

    ## compute weighted mean at last timestep
    if (nt == ntimes)
      coef(object,transform=transform) <- apply(params,1L,weighted.mean,w=weights)

    ## compute effective sample size, log-likelihood
    ## also do resampling if filtering has not failed
    xx <- try(
              .Call(
                    pfilter_computations,
                    x=X,
                    params=params,
                    Np=Np[nt+1],
                    rw=FALSE,
                    rw_sd=numeric(0),
                    predmean=FALSE,
                    predvar=FALSE,
                    filtmean=filter.mean,
                    onepar=FALSE,
                    weights=weights,
                    tol=tol
                    ),
              silent=FALSE
              )
    if (inherits(xx,'try-error')) {
      stop(sQuote("mif2.pfilter")," error",call.=FALSE)
    }
    all.fail <- xx$fail
    loglik[nt] <- xx$loglik
    eff.sample.size[nt] <- xx$ess

    x <- xx$states
    params <- xx$params
    if (filter.mean)
      filt.m[,nt] <- xx$fm

    if (all.fail) { ## all particles are lost
      nfail <- nfail+1
      if (verbose)
        message("filtering failure at time t = ",times[nt+1])
      if (nfail>max.fail)
        stop(sQuote("mif2.pfilter")," error: too many filtering failures",call.=FALSE)
    }

    if (verbose && (nt%%5==0))
      cat("mif2 pfilter timestep",nt,"of",ntimes,"finished\n")

  }

  if (nfail>0)
    warning(sprintf(ngettext(nfail,msg1="%d filtering failure occurred in ",
                             msg2="%d filtering failures occurred in "),nfail),
            sQuote("mif2.pfilter"),call.=FALSE)

  new(
      "pfilterd.pomp",
      object,
      paramMatrix=params,
      eff.sample.size=eff.sample.size,
      cond.loglik=loglik,
      filter.mean=filt.m,
      Np=Np,
      tol=tol,
      nfail=as.integer(nfail),
      loglik=sum(loglik)
      )
}

mif2.internal <- function (object, Nmif, start, Np, rw.sd, perturb.fn = NULL,
                           tol = 1e17, max.fail = Inf, transform = FALSE,
                           verbose = FALSE, .ndone = 0L,
                           .getnativesymbolinfo = TRUE, ...) {

  pompLoad(object)

  gnsi <- as.logical(.getnativesymbolinfo)

  Nmif <- as.integer(Nmif)
  if (Nmif<0) stop(sQuote("mif2")," error: ",sQuote("Nmif"),
                   " must be a positive integer",call.=FALSE)

  ntimes <- length(time(object))

  Np <- as.integer(Np)

  if (is.null(names(start)))
    stop(sQuote("mif2")," error: ",sQuote("start")," must be a named vector",
         call.=FALSE)

  cooling.fn <- cooling_fn(rw.sd,paramnames=names(start))

  conv.rec <- array(data=NA,dim=c(Nmif+1,length(start)+2),
                    dimnames=list(seq.int(.ndone,.ndone+Nmif),
                      c('loglik','nfail',names(start))))
  conv.rec[1L,] <- c(NA,NA,start)

  if (.ndone > 0) {                     # call is from 'continue'
    paramMatrix <- object@paramMatrix
  } else if (Nmif > 0) {                # initial call
    paramMatrix <- array(data=start,dim=c(length(start),Np[1L]),
                         dimnames=list(names(start),NULL))
  } else {                              # no work to do
    paramMatrix <- array(dim=c(0,0))
  }

  object <- as(object,"pomp")

  if (transform)
    paramMatrix <- partrans(object,paramMatrix,dir="toEstimationScale",
                            .getnativesymbolinfo=gnsi)

  ## iterate the filtering
  for (n in seq_len(Nmif)) {

    pfp <- try(
               mif2.pfilter(
                            object=object,
                            params=paramMatrix,
                            Np=Np,
                            mifiter=.ndone+n,
                            cooling.fn=cooling.fn,
                            perturb.fn=perturb.fn,
                            tol=tol,
                            max.fail=max.fail,
                            verbose=verbose,
                            filter.mean=(n==Nmif),
                            transform=transform,
                            .getnativesymbolinfo=gnsi
                            ),
               silent=FALSE
               )
    if (inherits(pfp,"try-error"))
      stop("mif2 particle-filter error")

    gnsi <- FALSE

    paramMatrix <- pfp@paramMatrix
    conv.rec[n+1,-c(1,2)] <- coef(pfp)
    conv.rec[n,c(1,2)] <- c(pfp@loglik,pfp@nfail)

    if (verbose) cat("mif2 iteration ",n," of ",Nmif," completed\n")

  }

  if (transform)
    pfp@paramMatrix <- partrans(object,paramMatrix,dir="fromEstimationScale",
                                .getnativesymbolinfo=gnsi)

  pompUnload(object)

  new(
      "mif2d.pomp",
      pfp,
      Nmif=Nmif,
      rw.sd=rw.sd,
      perturb.fn=perturb.fn,
      transform=transform,
      conv.rec=conv.rec,
      tol=tol
      )
}

setGeneric("mif2",function(object,...)standardGeneric("mif2"))

setMethod(
          "mif2",
          signature=signature(object="pomp"),
          definition = function (object, Nmif = 1, start, Np, rw.sd, perturb.fn,
            tol = 1e-17, max.fail = Inf, transform = FALSE,
            verbose = getOption("verbose"),...) {

            if (missing(start)) start <- coef(object)
            if (length(start)==0)
              stop(
                   sQuote("mif2")," error: ",sQuote("start")," must be specified if ",
                   sQuote("coef(object)")," is NULL",
                   call.=FALSE
                   )

            ntimes <- length(time(object))
            if (missing(Np)) {
              stop(sQuote("mif2")," error: ",sQuote("Np")," must be specified",call.=FALSE) }
            else if (is.function(Np)) {
              Np <- try(
                        vapply(seq.int(1,ntimes),Np,numeric(1)),
                        silent=FALSE
                        )
              if (inherits(Np,"try-error"))
                stop(sQuote("mif2")," error: if ",sQuote("Np"),
                     " is a function, it must return a single positive integer")
            } else if (!is.numeric(Np))
              stop(sQuote("mif2")," error: ",sQuote("Np"),
                   " must be a number, a vector of numbers, or a function")
            if (length(Np)==1) {
              Np <- rep(Np,times=ntimes+1)
            } else if (length(Np)==ntimes) {
              Np <- c(Np,Np[1L])
            } else if (length(Np)>ntimes) {
              if (Np[1L] != Np[ntimes+1])
                stop(sQuote("mif2")," error: Np[ntimes+1] != Np[1]")
              if (length(Np)>ntimes+1)
                warning("in ",sQuote("mif2"),": Np[k] ignored for k > ntimes")
            }
            if (any(Np<=0))
              stop("number of particles, ",sQuote("Np"),", must always be positive")

            if (missing(perturb.fn)) {
              perturb.fn <- function (theta, sd) {
                theta[names(sd),] <- rnorm(n=length(sd)*ncol(theta),mean=theta[names(sd),],sd=sd)
                theta
              }
            } else {
              perturb.fn <- match.fun(perturb.fn)
              if (!all(c('theta','sd')%in%names(formals(perturb.fn)))) {
                stop(
                     sQuote("mif2")," error: ",
                     sQuote("perturb.fn"),
                     " must be a function of prototype ",
                     sQuote("perturb.fn(theta,sd)"),
                     call.=FALSE
                     )
              }
            }

            mif2.internal(
                          object=object,
                          Nmif=Nmif,
                          start=start,
                          Np=Np,
                          rw.sd=rw.sd,
                          perturb.fn=perturb.fn,
                          tol=tol,
                          max.fail=max.fail,
                          transform=transform,
                          verbose=verbose,
                          ...
                          )

          }
          )


setMethod(
          "mif2",
          signature=signature(object="pfilterd.pomp"),
          definition = function (object, Nmif = 1, Np, tol, ...) {

            if (missing(Np)) Np <- object@Np
            if (missing(tol)) tol <- object@tol

            mif2(
                 object=as(object,"pomp"),
                 Nmif=Nmif,
                 Np=Np,
                 tol=tol,
                 ...
                 )
          }
          )

setMethod(
          "mif2",
          signature=signature(object="mif2d.pomp"),
          definition = function (object, Nmif, start, Np, rw.sd, perturb.fn, tol,
            transform, ...) {

            if (missing(Nmif)) Nmif <- object@Nmif
            if (missing(start)) start <- coef(object)
            if (missing(rw.sd)) rw.sd <- object@rw.sd
            if (missing(perturb.fn)) perturb.fn <- object@perturb.fn
            if (missing(transform)) transform <- object@transform

            if (missing(Np)) Np <- object@Np
            if (missing(tol)) tol <- object@tol

            f <- selectMethod("mif2","pomp")

            f(object,Nmif=Nmif,start=start,Np=Np,rw.sd=rw.sd,
              perturb.fn=perturb.fn,tol=tol,transform=transform,...)
          }
          )

setMethod(
          'continue',
          signature=signature(object='mif2d.pomp'),
          definition = function (object, Nmif = 1, ...) {

            ndone <- object@Nmif

            obj <- mif2(object=object,Nmif=Nmif,.ndone=ndone,...)

            object@conv.rec[ndone+1,c('loglik','nfail')] <- obj@conv.rec[1L,c('loglik','nfail')]
            obj@conv.rec <- rbind(
                                  object@conv.rec,
                                  obj@conv.rec[-1L,colnames(object@conv.rec)]
                                  )
            obj@Nmif <- as.integer(ndone+Nmif)

            obj
          }
          )

## extract the estimated log likelihood
setMethod('logLik','mif2d.pomp',function(object,...)object@loglik)

setMethod('conv.rec','mif2d.pomp',
          function (object, pars, transform = FALSE, ...) {
            conv.rec.internal(object=object,pars=pars,transform=transform,...)
          }
          )

## mif2List class
setClass(
         'mif2List',
         contains='list',
         validity=function (object) {
           if (!all(sapply(object,is,'mif2d.pomp'))) {
             retval <- paste0(
                              "error in ",sQuote("c"),
                              ": dissimilar objects cannot be combined"
                              )
             return(retval)
           }
           d <- sapply(object,function(x)dim(x@conv.rec))
           if (!all(apply(d,1,diff)==0)) {
             retval <- paste0(
                              "error in ",sQuote("c"),
                              ": to be combined, ",sQuote("mif2d.pomp"),
                              " objects must equal numbers of iterations"
                              )
             return(retval)
           }
           TRUE
         }
         )

setMethod(
          'c',
          signature=signature(x='mif2d.pomp'),
          definition=function (x, ...) {
            y <- list(...)
            if (length(y)==0) {
              new("mif2List",list(x))
            } else {
              p <- sapply(y,is,'mif2d.pomp')
              pl <- sapply(y,is,'mif2List')
              if (any(!(p||pl)))
                stop("cannot mix ",sQuote("mif2d.pomp"),
                     " and non-",sQuote("mif2d.pomp")," objects")
              y[p] <- lapply(y[p],list)
              y[pl] <- lapply(y[pl],as,"list")
              new("mif2List",c(list(x),y,recursive=TRUE))
            }
          }
          )

setMethod(
          'c',
          signature=signature(x='mif2List'),
          definition=function (x, ...) {
            y <- list(...)
            if (length(y)==0) {
              x
            } else {
              p <- sapply(y,is,'mif2d.pomp')
              pl <- sapply(y,is,'mif2List')
              if (any(!(p||pl)))
                stop("cannot mix ",sQuote("mif2d.pomp"),
                     " and non-",sQuote("mif2d.pomp")," objects")
              y[p] <- lapply(y[p],list)
              y[pl] <- lapply(y[pl],as,"list")
              new("mif2List",c(as(x,"list"),y,recursive=TRUE))
            }
          }
          )

setMethod(
          "[",
          signature=signature(x="mif2List"),
          definition=function(x, i, ...) {
            new('mif2List',as(x,"list")[i])
          }
          )

setMethod(
          'conv.rec',
          signature=signature(object='mif2List'),
          definition=function (object, ...) {
            lapply(object,conv.rec,...)
          }
          )

mif2.diagnostics <- function (z) {
  ## assumes that z is a list of mif2d.pomps with identical structure
  mar.multi <- c(0,5.1,0,2.1)
  oma.multi <- c(6,0,5,0)
  xx <- z[[1]]
  parnames <- names(coef(xx,transform=xx@transform))
  estnames <- parnames
  
  ## plot filter means
  filt.diag <- rbind("eff. sample size"=xx@eff.sample.size,filter.mean(xx))
  filtnames <- rownames(filt.diag)
  plotnames <- filtnames
  lognames <- filtnames[1] # eff. sample size
  nplots <- length(plotnames)
  n.per.page <- min(nplots,10)
  if(n.per.page<=4) nc <- 1 else nc <- 2
  nr <- ceiling(n.per.page/nc)
  oldpar <- par(mar=mar.multi,oma=oma.multi,mfcol=c(nr,nc),ask=dev.interactive(orNone=TRUE))
  on.exit(par(oldpar))
  low <- 1
  hi <- 0
  time <- time(xx)
  while (hi<nplots) {
    hi <- min(low+n.per.page-1,nplots)
    for (i in seq(from=low,to=hi,by=1)) {
      n <- i-low+1
      logplot <- if (plotnames[i]%in%lognames) "y" else ""
      dat <- sapply(
                    z,
                    function(po, label) {
                      if (label=="eff. sample size")
                        po@eff.sample.size
                      else
                        filter.mean(po,label)
                    },
                    label=plotnames[i]
                    )
      matplot(
              y=dat,
              x=time,
              axes = FALSE,
              xlab = "",
              log=logplot,
              ylab = "",
              type = "l"
              )
      box()
      y.side <- 2
      axis(y.side, xpd = NA)
      mtext(plotnames[i], y.side, line = 3)
      do.xax <- (n%%nr==0||n==n.per.page)
      if (do.xax) axis(1,xpd=NA)
      if (do.xax) mtext("time",side=1,line=3)
    }
    low <- hi+1
    mtext("Filter diagnostics (last iteration)",3,line=2,outer=TRUE)
  }

  ## plot mif convergence diagnostics
  other.diagnostics <- c("loglik", "nfail")
  plotnames <- c(other.diagnostics,estnames)
  nplots <- length(plotnames)
  n.per.page <- min(nplots,10)
  nc <- if (n.per.page<=4) 1 else 2
  nr <- ceiling(n.per.page/nc)
  par(mar=mar.multi,oma=oma.multi,mfcol=c(nr,nc))
  ## on.exit(par(oldpar))
  low <- 1
  hi <- 0
  iteration <- seq(0,xx@Nmif)
  while (hi<nplots) {
    hi <- min(low+n.per.page-1,nplots)
    for (i in seq(from=low,to=hi,by=1)) {
      n <- i-low+1
      dat <- sapply(z,function(po,label) conv.rec(po,label),label=plotnames[i])
      matplot(
              y=dat,
              x=iteration,
              axes = FALSE,
              xlab = "",
              ylab = "",
              type = "l"
              )
      box()
      y.side <- 2
      axis(y.side,xpd=NA)
      mtext(plotnames[i],y.side,line=3)
      do.xax <- (n%%nr==0||n==n.per.page)
      if (do.xax) axis(1,xpd=NA)
      if (do.xax) mtext("MIF iteration",side=1,line=3)
    }
    low <- hi+1
    mtext("MIF convergence diagnostics",3,line=2,outer=TRUE)
  }
  invisible(NULL)
}


setMethod(
          "plot",
          "mif2d.pomp",
          function (x, y, ...) {
            if (!missing(y)) {
              y <- substitute(y)
              warning(sQuote(y)," is ignored")
            }
            mif2.diagnostics(list(x))
          }
          )

setMethod(
          "plot",
          signature=signature(x='mif2List'),
          definition=function (x, y, ...) {
            if (!missing(y)) {
              y <- substitute(y)
              warning(sQuote(y)," is ignored")
            }
            mif2.diagnostics(x)
          }
          )
