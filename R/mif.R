## MIF algorithm functions

## define the mif class
setClass(
    'mif',
    contains='pfilterd.pomp',
    slots=c(
        transform = "logical",
        ivps = 'character',
        pars = 'character',
        Nmif = 'integer',
        particles = 'function',
        var.factor='numeric',
        ic.lag='integer',
        cooling.type='character',
        cooling.fraction.50='numeric',
        method='character',
        random.walk.sd = 'numeric',
        conv.rec = 'matrix'
    )
)

default.mif.particles.fun <- function (Np, center, sd, ...) {
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

mif.cooling.function <- function (type, perobs, fraction, ntimes) {
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
            if (fraction < 1) {
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
            } else {
                function (nt, m) {
                    list(alpha=1,gamma=1)
                }
            }
        },
        stop("unrecognized cooling schedule type ",sQuote(type),call.=FALSE)
    )
}

mif.pfilter <- function (object, params, Np,
                         tol, max.fail,
                         pred.mean = FALSE,
                         pred.var = FALSE,
                         cooling, cooling.m,
                         .mif2 = FALSE,
                         .rw.sd,
                         verbose = FALSE,
                         .transform = FALSE,
                         .getnativesymbolinfo = TRUE) {

    object <- as(object,"pomp")
    pompLoad(object)

    ptsi.for <- gnsi.rproc <- gnsi.dmeas <- as.logical(.getnativesymbolinfo)
    mif2 <- as.logical(.mif2)
    pred.mean <- as.logical(pred.mean)
    pred.var <- as.logical(pred.var)
    verbose <- as.logical(verbose)
    transform <- as.logical(.transform)

    one.par <- FALSE
    times <- time(object,t0=TRUE)
    ntimes <- length(times)-1

    if (missing(Np))
        Np <- NCOL(params)
    if (is.function(Np)) {
        Np <- try(
            vapply(seq.int(from=0,to=ntimes,by=1),Np,numeric(1)),
            silent=FALSE
        )
        if (inherits(Np,"try-error"))
            stop("if ",sQuote("Np")," is a function, it must return a single positive integer",call.=FALSE)
    }
    if (length(Np)==1)
        Np <- rep(Np,times=ntimes+1)
    else if (length(Np)!=(ntimes+1))
        stop(sQuote("Np")," must have length 1 or length ",ntimes+1,call.=FALSE)
    if (any(Np<=0))
        stop("number of particles, ",sQuote("Np"),", must always be positive",call.=FALSE)
    if (!is.numeric(Np))
        stop(sQuote("Np")," must be a number, a vector of numbers, or a function",call.=FALSE)
    Np <- as.integer(Np)

    if (is.null(dim(params))) {
        one.par <- TRUE               # there is only one parameter vector
        coef(object) <- params        # set params slot to the parameters
        params <- matrix(
            params,
            nrow=length(params),
            ncol=Np[1L],
            dimnames=list(
                names(params),
                NULL
            )
        )
    }
    paramnames <- rownames(params)
    if (is.null(paramnames))
        stop(sQuote("mif")," error: ",sQuote("params")," must have rownames",call.=FALSE)

    init.x <- init.state(
        object,
        params=if (transform) {
                   partrans(object,params,dir="fromEstimationScale",
                            .getnativesymbolinfo=ptsi.for)
               } else {
                   params
               }
    )
    statenames <- rownames(init.x)
    nvars <- nrow(init.x)
    ptsi.for <- FALSE
    x <- init.x

    random.walk <- !missing(.rw.sd)
    if (random.walk) {
        rw.names <- names(.rw.sd)
        if (is.null(rw.names)||!is.numeric(.rw.sd))
            stop(sQuote("mif")," error: ",sQuote(".rw.sd")," must be a named vector",call.=FALSE)
        if (!all(rw.names%in%paramnames))
            stop(
                sQuote("mif")," error: the rownames of ",
                sQuote("params")," must include all of the names of ",
                sQuote(".rw.sd"),"",call.=FALSE
            )
        sigma <- .rw.sd
    } else {
        rw.names <- character(0)
        sigma <- NULL
    }

    loglik <- rep(NA,ntimes)
    eff.sample.size <- numeric(ntimes)
    nfail <- 0
    npars <- length(rw.names)

    ## set up storage for prediction means, variances, etc.
    if (pred.mean)
        pred.m <- matrix(
        data=0,
        nrow=nvars+npars,
        ncol=ntimes,
        dimnames=list(
            variable=c(statenames,rw.names),
            time=time(object))
    )
    else
        pred.m <- array(data=numeric(0),dim=c(0,0))

    if (pred.var)
        pred.v <- matrix(
        data=0,
        nrow=nvars+npars,
        ncol=ntimes,
        dimnames=list(
            variable=c(statenames,rw.names),
            time=time(object))
    )
    else
        pred.v <- array(data=numeric(0),dim=c(0,0))

    if (random.walk) {
        filt.m <- matrix(
            data=0,
            nrow=nvars+length(paramnames),
            ncol=ntimes,
            dimnames=list(
                variable=c(statenames,paramnames),
                time=time(object)
            )
        )
    } else {
        filt.m <- matrix(
            data=0,
            nrow=nvars,
            ncol=ntimes,
            dimnames=list(
                variable=statenames,
                time=time(object))
        )
    }

    for (nt in seq_len(ntimes)) { ## main loop

        if (mif2) {
            cool.sched <- cooling(nt=nt,m=cooling.m)
            sigma1 <- sigma*cool.sched$alpha
        } else {
            sigma1 <- sigma
        }

        ## transform the parameters if necessary
        if (transform) tparams <- partrans(object,params,dir="fromEstimationScale",
                                           .getnativesymbolinfo=ptsi.for)
        ptsi.for <- FALSE

        ## advance the state variables according to the process model
        X <- try(
            rprocess(
                object,
                xstart=x,
                times=times[c(nt,nt+1)],
                params=if (transform) tparams else params,
                offset=1,
                .getnativesymbolinfo=gnsi.rproc
            ),
            silent=FALSE
        )
        if (inherits(X,'try-error'))
            stop(sQuote("mif")," error: process simulation error",call.=FALSE)
        gnsi.rproc <- FALSE

        if (pred.var) { ## check for nonfinite state variables and parameters
            problem.indices <- unique(which(!is.finite(X),arr.ind=TRUE)[,1L])
            if (length(problem.indices)>0) {  # state variables
                stop(
                    sQuote("mif")," error: non-finite state variable(s): ",
                    paste(rownames(X)[problem.indices],collapse=', '),
                    call.=FALSE
                )
            }
            if (random.walk) { # parameters (need to be checked only if 'random.walk=TRUE')
                problem.indices <- unique(which(!is.finite(params[rw.names,,drop=FALSE]),arr.ind=TRUE)[,1L])
                if (length(problem.indices)>0) {
                    stop(
                        sQuote("mif")," error: non-finite parameter(s): ",
                        paste(rw.names[problem.indices],collapse=', '),
                        call.=FALSE
                    )
                }
            }
        }

        ## determine the weights
        weights <- try(
            dmeasure(
                object,
                y=object@data[,nt,drop=FALSE],
                x=X,
                times=times[nt+1],
                params=if (transform) tparams else params,
                log=FALSE,
                .getnativesymbolinfo=gnsi.dmeas
            ),
            silent=FALSE
        )
        if (inherits(weights,'try-error'))
            stop("in ",sQuote("mif"),": error in calculation of weights.",call.=FALSE)
        if (!all(is.finite(weights)))
            stop("in ",sQuote("mif"),": ",sQuote("dmeasure")," returns non-finite value.",call.=FALSE)
        gnsi.dmeas <- FALSE

        ## compute prediction mean, prediction variance, filtering mean,
        ## effective sample size, log-likelihood
        ## also do resampling if filtering has not failed
        xx <- try(
            .Call(
                pfilter_computations1,
                x=X,
                params=params,
                Np=Np[nt+1],
                rw=random.walk,
                rw_sd=sigma1,
                predmean=pred.mean,
                predvar=pred.var,
                filtmean=TRUE,
                trackancestry=FALSE,
                onepar=one.par,
                weights=weights,
                tol=tol
            ),
            silent=FALSE
        )
        if (inherits(xx,'try-error')) {
            stop(sQuote("mif")," error",call.=FALSE)
        }
        all.fail <- xx$fail
        loglik[nt] <- xx$loglik
        eff.sample.size[nt] <- xx$ess

        x <- xx$states
        params <- xx$params

        if (pred.mean)
            pred.m[,nt] <- xx$pm
        if (pred.var)
            pred.v[,nt] <- xx$pv
        filt.m[,nt] <- xx$fm

        if (all.fail) { ## all particles are lost
            nfail <- nfail+1
            if (verbose)
                message("filtering failure at time t = ",times[nt+1])
            if (nfail>max.fail)
                stop(sQuote("mif")," error: too many filtering failures",call.=FALSE)
        }

        if (verbose && (nt%%5==0))
            cat("pfilter timestep",nt,"of",ntimes,"finished\n")

    } ## end of main loop

    if (nfail>0)
        warning(sprintf(ngettext(nfail,msg1="%d filtering failure occurred in ",
                                 msg2="%d filtering failures occurred in "),nfail),
                sQuote("mif"),call.=FALSE)

    pompUnload(object)

    new(
        "pfilterd.pomp",
        object,
        pred.mean=pred.m,
        pred.var=pred.v,
        filter.mean=filt.m,
        paramMatrix=if (mif2) params else array(data=numeric(0),dim=c(0,0)),
        eff.sample.size=eff.sample.size,
        cond.loglik=loglik,
        Np=as.integer(Np),
        tol=tol,
        nfail=as.integer(nfail),
        loglik=sum(loglik)
    )
}

mif.internal <- function (object, Nmif,
                          start, ivps,
                          particles,
                          rw.sd,
                          Np, var.factor, ic.lag,
                          cooling.type, cooling.fraction.50,
                          method,
                          tol, max.fail,
                          verbose, transform, .ndone = 0L,
                          paramMatrix = NULL,
                          .getnativesymbolinfo = TRUE,
                          ...) {

    if (method=="mif2") {
        warning(
            "method=",sQuote("mif2")," is deprecated and will be removed in a future release.\n",
            "Use ",sQuote("mif2")," instead.",call.=FALSE)
    }

    pompLoad(object)

    gnsi <- as.logical(.getnativesymbolinfo)

    transform <- as.logical(transform)

    if (length(start)==0)
        stop(
            "mif error: ",sQuote("start")," must be specified if ",
            sQuote("coef(object)")," is NULL",
            call.=FALSE
        )

    if (transform)
        start <- partrans(object,start,dir="toEstimationScale")

    start.names <- names(start)
    if (is.null(start.names))
        stop("mif error: ",sQuote("start")," must be a named vector",call.=FALSE)

    rw.names <- names(rw.sd)
    if (is.null(rw.names) || any(rw.sd<0))
        stop("mif error: ",sQuote("rw.sd")," must be a named non-negative numerical vector",call.=FALSE)
    if (!all(rw.names%in%start.names))
        stop("mif error: all the names of ",sQuote("rw.sd")," must be names of ",sQuote("start"),call.=FALSE)
    rw.names <- rw.names[rw.sd>0]
    rw.sd <- rw.sd[rw.sd>0]
    if (length(rw.names) == 0)
        stop("mif error: ",sQuote("rw.sd")," must have one positive entry for each parameter to be estimated",call.=FALSE)

    pars <- rw.names[!(rw.names%in%ivps)]
    
    if (!is.character(ivps) || !all(ivps%in%start.names))
        stop("mif error: ",sQuote("ivps")," must name model parameters",call.=FALSE)

    ntimes <- length(time(object))
    if (is.null(Np)) stop("mif error: ",sQuote("Np")," must be specified",call.=FALSE)
    if (is.function(Np)) {
        Np <- try(
            vapply(seq.int(from=0,to=ntimes,by=1),Np,numeric(1)),
            silent=FALSE
        )
        if (inherits(Np,"try-error"))
            stop("if ",sQuote("Np")," is a function, it must return a single positive integer",call.=FALSE)
    }
    if (length(Np)==1)
        Np <- rep(Np,times=ntimes+1)
    else if (length(Np)!=(ntimes+1))
        stop(sQuote("Np")," must have length 1 or length ",ntimes+1)
    if (any(Np<=0))
        stop("number of particles, ",sQuote("Np"),", must always be positive",call.=FALSE)
    if (!is.numeric(Np))
        stop(sQuote("Np")," must be a number, a vector of numbers, or a function",call.=FALSE)
    Np <- as.integer(Np)

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

    if (missing(cooling.fraction.50))
        stop("mif error: ",sQuote("cooling.fraction.50")," must be specified",call.=FALSE)
    cooling.fraction.50 <- as.numeric(cooling.fraction.50)
    if ((length(cooling.fraction.50)!=1)||(cooling.fraction.50<0)||(cooling.fraction.50>1))
        stop("mif error: ",sQuote("cooling.fraction.50")," must be a number between 0 and 1",call.=FALSE)

    cooling <- mif.cooling.function(
        type=cooling.type,
        perobs=(method=="mif2"),
        fraction=cooling.fraction.50,
        ntimes=ntimes
    )

    if ((method=="mif2")&&(Np[1L]!=Np[ntimes+1]))
        stop("the first and last values of ",sQuote("Np")," must agree when method = ",sQuote("mif2"),call.=FALSE)

    if ((length(var.factor)!=1)||(var.factor < 0))
        stop("mif error: ",sQuote("var.factor")," must be a positive number",call.=FALSE)

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
            iteration=seq(.ndone,.ndone+Nmif),
            variable=c('loglik','nfail',names(theta))
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

    have.parmat <- !(is.null(paramMatrix) || length(paramMatrix)==0)

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
            mif.pfilter(
                object=obj,
                params=P,
                Np=Np,
                tol=tol,
                max.fail=max.fail,
                pred.mean=(n==Nmif),
                pred.var=((method=="mif")||(n==Nmif)),
                cooling=cooling,
                cooling.m=.ndone+n,
                .mif2=(method=="mif2"),
                .rw.sd=sigma.n[pars],
                .transform=transform,
                verbose=verbose,
                .getnativesymbolinfo=gnsi
            ),
            silent=TRUE
        )
        if (inherits(pfp,"try-error"))
            stop("in ",sQuote("mif"),": error in ",sQuote("mif.pfilter"),
                 ":\n",pfp,call.=FALSE)

        gnsi <- FALSE

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
            stop("unrecognized method ",sQuote(method),call.=FALSE)
        )
        theta[ivps] <- pfp@filter.mean[ivps,ic.lag]
        conv.rec[n+1,-c(1,2)] <- theta
        conv.rec[n,c(1,2)] <- c(pfp@loglik,pfp@nfail)

        if (verbose) cat("MIF iteration ",n," of ",Nmif," completed\n")

    } ### end of main loop

    ## back transform the parameter estimate if necessary
    if (transform) theta <- partrans(pfp,theta,dir="fromEstimationScale")

    pompUnload(object)

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
        cooling.fraction.50=cooling.fraction.50,
        paramMatrix=if (method=="mif2") paramMatrix else array(data=numeric(0),dim=c(0,0))
    )
}

setMethod(
    "mif",
    signature=signature(object="pomp"),
    function (object, Nmif = 1,
              start,
              ivps = character(0),
              particles, rw.sd,
              Np, ic.lag, var.factor = 1,
              cooling.type = c("geometric","hyperbolic"),
              cooling.fraction.50,
              method = c("mif","unweighted","fp","mif2"),
              tol = 1e-17, max.fail = Inf,
              verbose = getOption("verbose"),
              transform = FALSE,
              ...) {

        method <- match.arg(method)

        if (missing(start)) start <- coef(object)
        if (missing(rw.sd))
            stop("mif error: ",sQuote("rw.sd")," must be specified",call.=FALSE)
        if (missing(ic.lag)) {
            if (length(ivps)>0 && (method != "mif2")) {
                stop("mif error: ",sQuote("ic.lag"),
                     " must be specified if ",sQuote("ivps"),
                     " are",call.=FALSE)
            } else {
                ic.lag <- length(time(object))
            }
        }

        if (missing(Np))
            stop("mif error: ",sQuote("Np")," must be specified",call.=FALSE)

        cooling.type <- match.arg(cooling.type)

        if (missing(particles)) { # use default: normal distribution
            particles <- default.mif.particles.fun
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
            ivps=ivps,
            particles=particles,
            rw.sd=rw.sd,
            Np=Np,
            cooling.type=cooling.type,
            cooling.fraction.50=cooling.fraction.50,
            var.factor=var.factor,
            ic.lag=ic.lag,
            method=method,
            tol=tol,
            max.fail=max.fail,
            verbose=verbose,
            transform=transform,
            ...
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
              ivps,
              particles, rw.sd,
              Np, ic.lag, var.factor,
              cooling.type, cooling.fraction.50,
              method,
              tol,
              transform,
              ...) {

        if (missing(Nmif)) Nmif <- object@Nmif
        if (missing(start)) start <- coef(object)
        if (missing(ivps)) ivps <- object@ivps
        if (missing(particles)) particles <- object@particles
        if (missing(rw.sd)) rw.sd <- object@random.walk.sd
        if (missing(ic.lag)) ic.lag <- object@ic.lag
        if (missing(var.factor)) var.factor <- object@var.factor
        if (missing(cooling.type)) cooling.type <- object@cooling.type
        if (missing(cooling.fraction.50)) cooling.fraction.50 <- object@cooling.fraction.50
        if (missing(method)) method <- object@method
        if (missing(transform)) transform <- object@transform

        if (missing(Np)) Np <- object@Np
        if (missing(tol)) tol <- object@tol

        mif(
            object=as(object,"pomp"),
            Nmif=Nmif,
            start=start,
            ivps=ivps,
            particles=particles,
            rw.sd=rw.sd,
            Np=Np,
            cooling.type=cooling.type,
            cooling.fraction.50=cooling.fraction.50,
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
        names(dimnames(obj@conv.rec)) <- c("iteration","variable")
        obj@Nmif <- as.integer(ndone+Nmif)

        obj
    }
)
