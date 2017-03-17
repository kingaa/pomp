## MIF2 algorithm functions

rw.sd <- safecall

pkern.sd <- function (rw.sd, time, paramnames) {
    ep <- paste0("in ",sQuote("mif2"),": ")
    if (is.matrix(rw.sd)) return(rw.sd)
    if (is(rw.sd,"safecall")) {
        enclos <- rw.sd@envir
        rw.sd <- as.list(rw.sd@call)[-1L]
    } else {
        stop(ep,sQuote("rw.sd")," should be specified using the ",sQuote("rw.sd"),
             " function. See ",sQuote("?mif2"),".",call.=FALSE)
    }
    if (is.null(names(rw.sd)) | any(names(rw.sd)==""))
        stop(ep,"in ",sQuote("rw.sd"),": parameters must be referenced by name.",call.=FALSE)
    if (!all(names(rw.sd) %in% paramnames)) {
        unrec <- names(rw.sd)[!names(rw.sd) %in% paramnames]
        stop(ep,"the following parameter(s), ",
             "given random walks in ",sQuote("rw.sd"),", are not present in ",
             sQuote("start"),": ",paste(sapply(unrec,sQuote),collapse=","),
             call.=FALSE)
    }
    ivp <- function (sd, lag = 1L) {
        sd*(seq_along(time)==lag)
    }
    sds <- lapply(rw.sd,eval,envir=list(time=time,ivp=ivp),enclos=enclos)
    for (n in names(sds)) {
        len <- length(sds[[n]])
        if (len==1) {
            sds[[n]] <- rep(sds[[n]],length(time))
        } else if (len!=length(time)) {
            stop(ep,sQuote("rw.sd")," spec for parameter ",sQuote(n),
                 " does not evaluate to a vector of the correct length (",
                 length(time),").",call.=FALSE)
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

mif2.cooling <- function (type, fraction, ntimes) {
    switch(
        type,
        geometric={
            factor <- fraction^(1/50)
            function (nt, m) {
                alpha <- factor^(nt/ntimes+m-1)
                list(alpha=alpha,gamma=alpha^2)
            }
        },
        hyperbolic={
            if (fraction < 1) {
                scal <- (50*ntimes*fraction-1)/(1-fraction)
                function (nt, m) {
                    alpha <- (1+scal)/(scal+nt+ntimes*(m-1))
                    list(alpha=alpha,gamma=alpha^2)
                }
            } else {
                function (nt, m) {
                    list(alpha=1,gamma=1)
                }
            }
        }
    )
}

mif2.pfilter <- function (object, params, Np,
                          mifiter, rw.sd, cooling.fn,
                          tol = 1e-17, max.fail = Inf,
                          transform, verbose,
                          .indices = integer(0),
                          .getnativesymbolinfo = TRUE) {

    ep <- paste0("in ",sQuote("mif2.pfilter"),": ")

    gnsi <- as.logical(.getnativesymbolinfo)
    transform <- as.logical(transform)
    verbose <- as.logical(verbose)
    mifiter <- as.integer(mifiter)
    Np <- as.integer(Np)

    do_ta <- length(.indices)>0L
    if (do_ta && length(.indices)!=Np[1L])
        stop(ep,sQuote(".indices"),
             " has improper length",call.=FALSE)

    times <- time(object,t0=TRUE)
    ntimes <- length(times)-1

    loglik <- rep(NA,ntimes)
    eff.sample.size <- numeric(ntimes)
    nfail <- 0

    for (nt in seq_len(ntimes)) {

        ## perturb parameters
        pmag <- cooling.fn(nt,mifiter)$alpha*rw.sd[,nt]
        params <- .Call(randwalk_perturbation,params,pmag)

        if (transform)
            tparams <- partrans(object,params,dir="fromEstimationScale",
                                .getnativesymbolinfo=gnsi)

        if (nt == 1L) {
            ## get initial states
            x <- init.state(object,params=if (transform) tparams else params)
        }

        ## advance the state variables according to the process model
        X <- tryCatch(
            rprocess(
                object,
                xstart=x,
                times=times[c(nt,nt+1)],
                params=if (transform) tparams else params,
                offset=1,
                .getnativesymbolinfo=gnsi
            ),
            error = function (e) {
                stop(ep,"process simulation error: ",
                     conditionMessage(e),call.=FALSE)
            }
        )

        ## determine the weights
        weights <- tryCatch(
            dmeasure(
                object,
                y=object@data[,nt,drop=FALSE],
                x=X,
                times=times[nt+1],
                params=if (transform) tparams else params,
                log=FALSE,
                .getnativesymbolinfo=gnsi
            ),
            error = function (e) {
                stop(ep,"error in calculation of weights: ",
                     conditionMessage(e),call.=FALSE)
            }
        )
        if (!all(is.finite(weights))) {
            first <- which(!is.finite(weights))[1L]
            datvals <- object@data[,nt]
            weight <- weights[first]
            states <- X[,first,1L]
            params <- if (transform) tparams[,first] else params[,first]
            msg <- paste0(
                sQuote("dmeasure")," returns non-finite value.\n",
                "likelihood, data, states, and parameters are:\n",
                "time: ",times[nt+1],"\n",
                "lik: ",weight,"\n",
                paste0(names(datvals),": ",datvals,collapse="\n"),"\n",
                paste0(names(states),": ",states,collapse="\n"),"\n",
                paste0(names(params),": ",params,collapse="\n")
            )
            stop(ep,msg,call.=FALSE)
        }
        gnsi <- FALSE

        ## compute weighted mean at last timestep
        if (nt == ntimes) {
            if (any(weights>0)) {
                coef(object,transform=transform) <- apply(params,1L,weighted.mean,w=weights)
            } else {
                warning(ep,"filtering failure at last filter iteration, using unweighted mean for ",
                        sQuote("coef"),call.=FALSE)
                coef(object,transform=transform) <- apply(params,1L,mean)
            }
        }

        ## compute effective sample size, log-likelihood
        ## also do resampling if filtering has not failed
        xx <- tryCatch(
            .Call(
                pfilter_computations,
                x=X,
                params=params,
                Np=Np[nt+1],
                rw_sd=numeric(0),
                predmean=FALSE,
                predvar=FALSE,
                filtmean=FALSE,
                trackancestry=do_ta,
                onepar=FALSE,
                weights=weights,
                tol=tol
            ),
            error = function (e) {
                stop(ep,"particle-filter error: ",conditionMessage(e),call.=FALSE) # nocov
            }
        )
        all.fail <- xx$fail
        loglik[nt] <- xx$loglik
        eff.sample.size[nt] <- xx$ess
        if (do_ta) {
            .indices <- .indices[xx$ancestry]
        }

        x <- xx$states
        params <- xx$params

        if (all.fail) { ## all particles are lost
            nfail <- nfail+1
            if (verbose)
                message("filtering failure at time t = ",times[nt+1])
            if (nfail>max.fail)
                stop(ep,"too many filtering failures",call.=FALSE)
        }

        if (verbose && (nt%%5==0))
            cat("mif2 pfilter timestep",nt,"of",ntimes,"finished\n")

    }

    if (nfail>0) {
        warning(
            ep,nfail,
            ngettext(
                nfail,
                msg1=" filtering failure occurred.",
                msg2=" filtering failures occurred."
            ),
            call.=FALSE
        )
    }

    new(
        "pfilterd.pomp",
        object,
        paramMatrix=params,
        eff.sample.size=eff.sample.size,
        cond.loglik=loglik,
        indices=.indices,
        Np=Np,
        tol=tol,
        nfail=as.integer(nfail),
        loglik=sum(loglik)
    )
}

mif2.internal <- function (object, Nmif, start, Np, rw.sd, transform = FALSE,
                           cooling.type, cooling.fraction.50,
                           tol = 1e-17, max.fail = Inf,
                           verbose = FALSE, .ndone = 0L,
                           .indices = integer(0),
                           .paramMatrix = NULL,
                           .getnativesymbolinfo = TRUE, ...) {

    ep <- paste0("in ",sQuote("mif2"),": ")

    pompLoad(object,verbose=verbose)

    transform <- as.logical(transform)
    verbose <- as.logical(verbose)
    gnsi <- as.logical(.getnativesymbolinfo)
    Np <- c(Np,Np[1L])

    if (Nmif <= 0)
        stop(ep,sQuote("Nmif")," must be a positive integer",call.=FALSE)

    cooling.fn <- mif2.cooling(
        type=cooling.type,
        fraction=cooling.fraction.50,
        ntimes=length(time(object))
    )

    if (is.null(.paramMatrix)) {
        if (.ndone > 0) {               # call is from 'continue'
            paramMatrix <- object@paramMatrix
            start <- apply(paramMatrix,1L,mean)
        } else {                         # initial call
            paramMatrix <- array(data=start,dim=c(length(start),Np[1L]),
                                 dimnames=list(variable=names(start),rep=NULL))
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

        pfp <- tryCatch(
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
                transform=transform,
                .indices=.indices,
                .getnativesymbolinfo=gnsi
            ),
            error = function (e) {
                stop(ep,conditionMessage(e),call.=FALSE)
            }
        )

        gnsi <- FALSE

        paramMatrix <- pfp@paramMatrix
        conv.rec[n+1,-c(1,2)] <- coef(pfp)
        conv.rec[n,c(1,2)] <- c(pfp@loglik,pfp@nfail)
        .indices <- pfp@indices

        if (verbose) cat("mif2 iteration",n,"of",Nmif,"completed\n")

    }

    if (transform)
        pfp@paramMatrix <- partrans(object,paramMatrix,dir="fromEstimationScale",
                                    .getnativesymbolinfo=gnsi)

    pompUnload(object,verbose=verbose)

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

        ep <- paste0("in ",sQuote("mif2"),": ")

        Nmif <- as.integer(Nmif)

        if (missing(start)) start <- coef(object)
        if (length(start)==0)
            stop(ep,sQuote("start")," must be specified if ",
                 sQuote("coef(object)")," is NULL",call.=FALSE)
        if (is.null(names(start)))
            stop(ep,sQuote("start")," must be a named vector",
                 call.=FALSE)

        ntimes <- length(time(object))

        if (missing(Np))
            stop(ep,sQuote("Np")," must be specified",call.=FALSE)
        else if (is.function(Np)) {
            Np <- tryCatch(
                vapply(seq_len(ntimes),Np,numeric(1)),
                error = function (e) {
                    stop(ep,"if ",sQuote("Np"),
                         " is a function, it must return a single positive integer",
                         call.=FALSE)
                }
            )
        } else if (!is.numeric(Np))
            stop(ep,sQuote("Np"),
                 " must be a number, a vector of numbers, or a function",
                 call.=FALSE)
        if (length(Np)==1) {
            Np <- rep(Np,times=ntimes)
        } else if (length(Np)>ntimes) {
            if (Np[1L] != Np[ntimes+1] || length(Np) > ntimes+1) {
                warning(ep,"Np[k] ignored for k > ntimes",call.=FALSE)
            }
            Np <- head(Np,ntimes)
        }
        if (any(Np <= 0))
            stop(ep,"number of particles, ",
                 sQuote("Np"),", must always be positive",call.=FALSE)

        if (missing(rw.sd))
            stop(ep,sQuote("rw.sd")," must be specified!",call.=FALSE)
        rw.sd <- pkern.sd(rw.sd,time=time(object),paramnames=names(start))

        cooling.type <- match.arg(cooling.type)

        cooling.fraction.50 <- as.numeric(cooling.fraction.50)
        if (cooling.fraction.50 <= 0 || cooling.fraction.50 > 1)
            stop(ep,sQuote("cooling.fraction.50"),
                 " must be in (0,1]",call.=FALSE)

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

        f <- selectMethod("mif2","mif2d.pomp")
        obj <- f(object=object,Nmif=Nmif,.ndone=ndone,...)

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
