## MIF2 algorithm functions

rw.sd <- function (...) {
    as.list(match.call())[-1L]
}

pkern.sd <- function (rw.sd, time, paramnames) {
    if (!all(names(rw.sd) %in% paramnames)) {
        unrec <- names(rw.sd)[!names(rw.sd) %in% paramnames]
        stop(sQuote("mif2")," error: the following parameter(s), ",
             "which are supposed to be estimated, are not present: ",
             paste(sapply(sQuote,unrec),collapse=","),
             call.=FALSE)
    }
    ivp <- function (sd, lag = 1L) {
        sd*(seq_along(time)==lag)
    }
    sds <- lapply(rw.sd,eval,envir=list(time=time,ivp=ivp))
    for (n in names(sds)) {
        len <- length(sds[[n]])
        if (len==1) {
            sds[[n]] <- rep(sds[[n]],length(time))
        } else if (len!=length(time)) {
            stop(sQuote("mif2")," error: ",sQuote("rw.sd"),
                 " spec for parameter ",sQuote(n),
                 " does not evaluate to a vector of the correct length (",
                 length(time),")",call.=FALSE)
        }
    }
    do.call(rbind,sds)
}

## define the mif2d.pomp class
setClass(
    'mif2d.pomp',
    contains='pfilterd.pomp',
    slots=c(
        Nmif = 'integer',
        rw.sd = 'matrix',
        cooling.type = 'character',
        cooling.fraction.50 = 'numeric',
        transform = 'logical',
        conv.rec = 'matrix'
    )
)

mif2.pfilter <- function (object, params, Np,
                          mifiter, rw.sd, cooling.fn,
                          tol = 1e-17, max.fail = Inf,
                          transform, verbose, filter.mean,
                          .getnativesymbolinfo = TRUE) {

    gnsi <- as.logical(.getnativesymbolinfo)
    transform <- as.logical(transform)
    verbose <- as.logical(verbose)
    filter.mean <- as.logical(filter.mean)
    mifiter <- as.integer(mifiter)
    Np <- as.integer(Np)

    times <- time(object,t0=TRUE)
    ntimes <- length(times)-1

    loglik <- rep(NA,ntimes)
    eff.sample.size <- numeric(ntimes)
    nfail <- 0

    for (nt in seq_len(ntimes)) {

        ## perturb parameters
        pmag <- cooling.fn(nt,mifiter)$alpha*rw.sd[,nt]
        params <- .Call(mif2_computations,params,pmag)

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
            silent=TRUE
        )
        if (inherits(X,'try-error'))
            stop("in ",sQuote("mif2.pfilter"),": process simulation error:",
                 X,call.=FALSE)

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
            silent=TRUE
        )
        if (inherits(weights,'try-error'))
            stop("in ",sQuote("mif2.pfilter"),": error in calculation of weights: ",
                 weights,call.=FALSE)
        if (!all(is.finite(weights)))
            stop(sQuote("mif2.pfilter")," error: ",sQuote("dmeasure"),
                 " returns non-finite value",call.=FALSE)
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
                trackancestry=FALSE,
                onepar=FALSE,
                weights=weights,
                tol=tol
            ),
            silent=TRUE
        )
        if (inherits(xx,'try-error')) {
            stop("in ",sQuote("mif2.pfilter"),": pfilter computation error: ",
                 xx,call.=FALSE)
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

mif2.internal <- function (object, Nmif, start, Np, rw.sd, transform = FALSE,
                           cooling.type, cooling.fraction.50,
                           tol = 1e17, max.fail = Inf, 
                           verbose = FALSE, .ndone = 0L,
                           .paramMatrix = NULL, 
                           .getnativesymbolinfo = TRUE, ...) {
    
    pompLoad(object)
    
    transform <- as.logical(transform)
    verbose <- as.logical(verbose)
    gnsi <- as.logical(.getnativesymbolinfo)
    Np <- c(Np,Np[1L])

    cooling.fn <- mif.cooling.function(
        type=cooling.type,
        perobs=TRUE,
        fraction=cooling.fraction.50,
        ntimes=length(time(object))
    )

    if (is.null(.paramMatrix)) {
        if (.ndone > 0) {              # call is from 'continue'
            paramMatrix <- object@paramMatrix
            start <- apply(paramMatrix,1L,mean)
        } else if (Nmif > 0) {         # initial call
            paramMatrix <- array(data=start,dim=c(length(start),Np[1L]),
                                 dimnames=list(names(start),NULL))
        } else {                        # no work to do
            paramMatrix <- array(dim=c(0,0))
        }
    } else { 
        paramMatrix <- .paramMatrix
        start <- apply(paramMatrix,1L,mean)
    }

    conv.rec <- array(dim=c(Nmif+1,length(start)+2),
                      dimnames=list(iteration=seq.int(.ndone,.ndone+Nmif),
                                    variable=c('loglik','nfail',names(start))))
    conv.rec[1L,] <- c(NA,NA,start)
    
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
                rw.sd=rw.sd,
                tol=tol,
                max.fail=max.fail,
                verbose=verbose,
                filter.mean=(n==Nmif),
                transform=transform,
                .getnativesymbolinfo=gnsi
            ),
            silent=TRUE
        )
        if (inherits(pfp,"try-error"))
            stop("in ",sQuote("mif2"),": particle-filter error:",pfp,call.=FALSE)

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
        cooling.type=cooling.type,
        cooling.fraction.50=cooling.fraction.50,
        transform=transform,
        conv.rec=conv.rec
    )
}

setMethod(
    "mif2",
    signature=signature(object="pomp"),
    definition = function (object, Nmif = 1, start, Np,
                           rw.sd, transform = FALSE,
                           cooling.type = c("hyperbolic", "geometric"),
                           cooling.fraction.50,
                           tol = 1e-17, max.fail = Inf,
                           verbose = getOption("verbose"),...) {

        Nmif <- as.integer(Nmif)
        if (Nmif<0) stop(sQuote("mif2")," error: ",sQuote("Nmif"),
                         " must be a positive integer",call.=FALSE)

        if (missing(start)) start <- coef(object)
        if (length(start)==0)
            stop(
                sQuote("mif2")," error: ",sQuote("start")," must be specified if ",
                sQuote("coef(object)")," is NULL",
                call.=FALSE
            )
        if (is.null(names(start)))
            stop(sQuote("mif2")," error: ",sQuote("start")," must be a named vector",
                 call.=FALSE)

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
            Np <- rep(Np,times=ntimes)
        } else if (length(Np)>ntimes) {
            if (Np[1L] != Np[ntimes+1] || length(Np) > ntimes+1) {
                warning("in ",sQuote("mif2"),": Np[k] ignored for k > ntimes")
            }
            Np <- head(Np,ntimes)
        }
        if (any(Np <= 0))
            stop("number of particles, ",sQuote("Np"),", must always be positive")

        if (missing(rw.sd))
            stop(sQuote("mif2")," error: ",sQuote("rw.sd")," must be specified!",call.=FALSE)
        if (!is.matrix(rw.sd)) {
            rw.sd <- pkern.sd(rw.sd,time=time(object),paramnames=names(start))
        }

        cooling.type <- match.arg(cooling.type)

        cooling.fraction.50 <- as.numeric(cooling.fraction.50)
        if (cooling.fraction.50 <= 0 || cooling.fraction.50 > 1)
            stop(sQuote("mif2")," error: ",
                 sQuote("cooling.fraction.50")," must be in (0,1]",call.=FALSE)

        mif2.internal(
            object=object,
            Nmif=Nmif,
            start=start,
            Np=Np,
            rw.sd=rw.sd,
            transform=transform,
            cooling.type=cooling.type,
            cooling.fraction.50=cooling.fraction.50,
            tol=tol,
            max.fail=max.fail,
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

        f <- selectMethod("mif2","pomp")
        f(object=object,Nmif=Nmif,Np=Np,tol=tol,...)
    }
)

setMethod(
    "mif2",
    signature=signature(object="mif2d.pomp"),
    definition = function (object, Nmif, start, Np,
                           rw.sd, transform, cooling.type, cooling.fraction.50,
                           tol, ...) {

        if (missing(Nmif)) Nmif <- object@Nmif
        if (missing(start)) start <- coef(object)
        if (missing(rw.sd)) rw.sd <- object@rw.sd
        if (missing(transform)) transform <- object@transform
        if (missing(cooling.type)) cooling.type <- object@cooling.type
        if (missing(cooling.fraction.50)) cooling.fraction.50 <- object@cooling.fraction.50
        
        if (missing(Np)) Np <- object@Np
        if (missing(tol)) tol <- object@tol

        f <- selectMethod("mif2","pomp")

        f(object,Nmif=Nmif,start=start,Np=Np,rw.sd=rw.sd,transform=transform,
          cooling.type=cooling.type,cooling.fraction.50=cooling.fraction.50,
          tol=tol,...)
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
        names(dimnames(obj@conv.rec)) <- c("iteration","variable")
        obj@Nmif <- as.integer(ndone+Nmif)

        obj
    }
)
