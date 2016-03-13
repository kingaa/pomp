## basic constructor of the pomp class
pomp.constructor <- function (data, times, t0, rprocess, dprocess,
                              rmeasure, dmeasure, measurement.model,
                              skeleton,
                              skeleton.type = c("map","vectorfield"),
                              skelmap.delta.t,
                              initializer, rprior, dprior,
                              params, covar, tcovar,
                              obsnames, statenames, paramnames, covarnames,
                              zeronames, PACKAGE,
                              fromEstimationScale, toEstimationScale,
                              globals, cdir, cfile, userdata,
                              ...,
                              .solibs = list(),
                              verbose = getOption("verbose",FALSE)) {

    ## preliminary error checking
    if (missing(data)) stop(sQuote("data")," is a required argument")
    if (missing(times)) stop(sQuote("times")," is a required argument")
    if (missing(t0)) stop(sQuote("t0")," is a required argument")
    if (missing(params)) params <- numeric(0)

    if (missing(cdir)) cdir <- NULL
    if (missing(cfile)) cfile <- NULL
    if (missing(userdata)) userdata <- list()
    added.userdata <- list(...)
    if (length(added.userdata)>0) {
        message("In ",sQuote("pomp"),
                ": the following unrecognized argument(s) ",
                "will be stored for use by user-defined functions: ",
                paste(sQuote(names(added.userdata)),collapse=","))
        userdata[names(added.userdata)] <- added.userdata
    }

    ## name of shared object library
    if (missing(PACKAGE)) PACKAGE <- NULL
    PACKAGE <- as.character(PACKAGE)

    if (missing(globals)) globals <- NULL
    if (!is(globals,"Csnippet"))
        globals <- as.character(globals)

    ## deal with missing components
    if (missing(skeleton)) skeleton <- NULL
    if (missing(rmeasure)) rmeasure <- NULL
    if (missing(dmeasure)) dmeasure <- NULL
    if (missing(rprior)) rprior <- NULL
    if (missing(dprior)) dprior <- NULL
    if (missing(fromEstimationScale)) fromEstimationScale <- NULL
    if (missing(toEstimationScale)) toEstimationScale <- NULL

    ## defaults for names of states, parameters, and accumulator variables
    if (missing(statenames)) statenames <- character(0)
    if (missing(paramnames)) paramnames <- character(0)
    if (missing(zeronames)) zeronames <- character(0)

    ## check the parameters and force them to be double-precision
    if (length(params)>0) {
        if (is.null(names(params)) || !is.numeric(params))
            stop("pomp error: ",sQuote("params")," must be a named numeric vector")
    }
    storage.mode(params) <- 'double'

    ## check the data and store it as double-precision matrix
    if (!is.numeric(data))
        stop("pomp error: ",sQuote("data")," must be numeric")
    storage.mode(data) <- 'double'
    if (missing(obsnames) || length(obsnames)==0) obsnames <- rownames(data)
    obsnames <- as.character(obsnames)

    ## check times
    if (!is.numeric(times) || !all(diff(times)>0))
        stop("pomp error: ",sQuote("times")," must be an increasing numeric vector")
    storage.mode(times) <- 'double'

    ## check t0
    if (!is.numeric(t0) || length(t0) > 1)
        stop("pomp error: the zero-time ",sQuote("t0")," must be a single number")
    storage.mode(t0) <- 'double'

    ## check and arrange covariates
    if (missing(covar)) {
        covar <- matrix(data=0,nrow=0,ncol=0)
        tcovar <- numeric(0)
    } else if (missing(tcovar)) {
        stop("pomp error: if ",sQuote("covar")," is supplied, ",
             sQuote("tcovar")," must also be supplied")
    } else if (is.data.frame(covar)) {
        if ((is.numeric(tcovar) && (tcovar<1 || tcovar>length(covar))) ||
            (is.character(tcovar) && (!(tcovar%in%names(covar)))) ||
            (!is.numeric(tcovar) && !is.character(tcovar))) {
            stop("pomp error: if ",sQuote("covar")," is a data frame, ",
                 sQuote("tcovar")," should indicate the time variable")
        } else if (is.numeric(tcovar)) {
            tpos <- tcovar
            tcovar <- covar[[tpos]]
            covar <- as.matrix(covar[-tpos])
        } else if (is.character(tcovar)) {
            tpos <- match(tcovar,names(covar))
            tcovar <- covar[[tpos]]
            covar <- as.matrix(covar[-tpos])
        }
    } else {
        covar <- as.matrix(covar)
    }
    if (missing(covarnames) || length(covarnames)==0) covarnames <- as.character(colnames(covar))
    if (!all(covarnames%in%colnames(covar))) {
        missing <- covarnames[!(covarnames%in%colnames(covar))]
        stop("pomp error: covariate(s) ",paste(missing,collapse=","),
             " are not among the columns of ",sQuote("covar"))
    }
    storage.mode(tcovar) <- "double"
    storage.mode(covar) <- "double"

    ## handle initializer
    default.init <- missing(initializer) || (is(initializer,"pomp.fun") && initializer@mode == pompfunmode$undef)
    if (default.init) initializer <- pomp.fun()

    ## default rprocess & dprocess
    if (missing(rprocess))
        rprocess <- function (xstart,times,params,...) stop(sQuote("rprocess")," not specified")
    if (missing(dprocess))
        dprocess <- function (x,times,params,log=FALSE,...) stop(sQuote("dprocess")," not specified")

    ## handle C snippets
    snips <- list()
    if (is(rprocess,"pompPlugin") && rprocess@csnippet)
        snips <- c(snips,setNames(list(slot(rprocess,rprocess@slotname)@text),rprocess@slotname))
    if (is(dprocess,"pompPlugin") && dprocess@csnippet)
        snips <- c(snips,setNames(list(slot(dprocess,dprocess@slotname)@text),dprocess@slotname))
    if (is(skeleton,"Csnippet"))
        snips <- c(snips,skeleton=skeleton@text)
    if (is(rmeasure,"Csnippet"))
        snips <- c(snips,rmeasure=rmeasure@text)
    if (is(dmeasure,"Csnippet"))
        snips <- c(snips,dmeasure=dmeasure@text)
    if (is(rprior,"Csnippet"))
        snips <- c(snips,rprior=rprior@text)
    if (is(dprior,"Csnippet"))
        snips <- c(snips,dprior=dprior@text)
    if (is(fromEstimationScale,"Csnippet"))
        snips <- c(snips,fromEstimationScale=fromEstimationScale@text)
    if (is(toEstimationScale,"Csnippet"))
        snips <- c(snips,toEstimationScale=toEstimationScale@text)
    if (is(initializer,"Csnippet"))
        snips <- c(snips,initializer=initializer@text)
    if (length(snips)>0) {
        libname <- try(
            do.call(
                pompCBuilder,
                c(
                    list(
                        dir=cdir,
                        name=cfile,
                        obsnames=obsnames,
                        statenames=statenames,
                        paramnames=paramnames,
                        covarnames=covarnames,
                        globals=globals,
                        verbose=verbose
                    ),
                    snips
                )
            ),
            silent=FALSE
        )
        if (inherits(libname,"try-error")) {
            stop("in ",sQuote("pomp"),": error in building shared-object library from Csnippets:\n",
                 libname,call.=FALSE)
        } else {
            .solibs <- c(.solibs,list(libname))
            libname <- libname$name
        }
    } else {
        libname <- ''
    }

    ## handle initializer
    if (!default.init) {
        initializer <- pomp.fun(
            f=initializer,
            PACKAGE=PACKAGE,
            proto=quote(initializer(params,t0,...)),
            slotname="initializer",
            libname=libname,
            statenames=statenames,
            paramnames=paramnames,
            obsnames=obsnames,
            covarnames=covarnames
        )
    }

    ## handle rprocess
    if (is(rprocess,"pompPlugin")) {
        rprocess <- plugin.handler(
            rprocess,
            libname=libname,
            statenames=statenames,
            paramnames=paramnames,
            obsnames=obsnames,
            covarnames=covarnames
        )
    }
    if (!is.function(rprocess))
        stop("pomp error: ",sQuote("rprocess")," must be a function")

    ## handle dprocess
    if (is(dprocess,"pompPlugin")) {
        dprocess <- plugin.handler(
            dprocess,
            libname=libname,
            statenames=statenames,
            paramnames=paramnames,
            obsnames=obsnames,
            covarnames=covarnames
        )
    }
    if (!is.function(dprocess))
        stop("pomp error: ",sQuote("dprocess")," must be a function")

    ## handle skeleton
    skeleton <- pomp.fun(
        f=skeleton,
        PACKAGE=PACKAGE,
        proto=quote(skeleton(x,t,params,...)),
        slotname="skeleton",
        libname=libname,
        statenames=statenames,
        paramnames=paramnames,
        obsnames=obsnames,
        covarnames=covarnames
    )

    ## type of skeleton (map or vectorfield)
    ## skelmap.delta.t has no meaning in the latter case
    skeleton.type <- match.arg(skeleton.type)
    skelmap.delta.t <- as.numeric(skelmap.delta.t)
    if (skelmap.delta.t <= 0)
        stop(sQuote("skelmap.delta.t")," must be positive")

    ## if 'measurement model' is specified as a formula, this overrides
    ## specification of 'rmeasure' or 'dmeasure'
    if (!missing(measurement.model)) {
        if (!(is.null(dmeasure)&&is.null(rmeasure))) {
            warning(
                "specifying ",sQuote("measurement.model"),
                " overrides specification of ",
                sQuote("rmeasure")," and ",sQuote("dmeasure")
            )
        }
        mm <- measform2pomp(measurement.model)
        rmeasure <- mm$rmeasure
        dmeasure <- mm$dmeasure
    }

    ## handle rmeasure
    rmeasure <- pomp.fun(
        f=rmeasure,
        PACKAGE=PACKAGE,
        proto=quote(rmeasure(x,t,params,...)),
        slotname="rmeasure",
        libname=libname,
        statenames=statenames,
        paramnames=paramnames,
        obsnames=obsnames,
        covarnames=covarnames
    )

    ## handle dmeasure
    dmeasure <- pomp.fun(
        f=dmeasure,
        PACKAGE=PACKAGE,
        proto=quote(dmeasure(y,x,t,params,log,...)),
        slotname="dmeasure",
        libname=libname,
        statenames=statenames,
        paramnames=paramnames,
        obsnames=obsnames,
        covarnames=covarnames
    )

    ## handle rprior
    rprior <- pomp.fun(
        f=rprior,
        PACKAGE=PACKAGE,
        proto=quote(rprior(params,...)),
        slotname="rprior",
        libname=libname,
        statenames=statenames,
        paramnames=paramnames,
        obsnames=obsnames,
        covarnames=covarnames
    )

    ## handle dprior
    if (is.null(dprior)) {
        ## by default, use flat improper prior
        dprior <- pomp.fun(f="_pomp_default_dprior",PACKAGE="pomp")
    } else {
        dprior <- pomp.fun(
            f=dprior,
            PACKAGE=PACKAGE,
            proto=quote(dprior(params,log,...)),
            slotname="dprior",
            libname=libname,
            statenames=statenames,
            paramnames=paramnames,
            obsnames=obsnames,
            covarnames=covarnames
        )
    }

    ## handle parameter transformations
    mpt <- is.null(fromEstimationScale)
    mpit <- is.null(toEstimationScale)
    if (xor(mpt,mpit)) {
        stop("if one of ",sQuote("fromEstimationScale"),", ",
             sQuote("toEstimationScale"),
             " is supplied, then so must the other")
    }
    has.trans <- !mpt
    if (has.trans) {
        from.trans <- pomp.fun(
            f=fromEstimationScale,
            PACKAGE=PACKAGE,
            proto=quote(from.trans(params,...)),
            slotname="fromEstimationScale",
            libname=libname,
            statenames=statenames,
            paramnames=paramnames,
            obsnames=obsnames,
            covarnames=covarnames
        )
        to.trans <- pomp.fun(
            f=toEstimationScale,
            PACKAGE=PACKAGE,
            proto=quote(to.trans(params,...)),
            slotname="toEstimationScale",
            libname=libname,
            statenames=statenames,
            paramnames=paramnames,
            obsnames=obsnames,
            covarnames=covarnames
        )
    } else {
        from.trans <- pomp.fun()
        to.trans <- pomp.fun()
    }
    if (has.trans &&
        from.trans@mode==pompfunmode$undef &&
        to.trans@mode==pompfunmode$undef
        ) has.trans <- FALSE

    if (nrow(covar)>0) {
        if (
        (skeleton@mode==pompfunmode$Rfun)
        &&!("covars"%in%names(formals(skeleton@R.fun)))
        )
            warning("a covariate table has been given, yet the ",
                    sQuote("skeleton")," function does not have ",
                    sQuote("covars")," as a formal argument: see ",
                    sQuote("?pomp"),call.=FALSE)
        if (
        (rmeasure@mode==pompfunmode$Rfun)
        &&!("covars"%in%names(formals(rmeasure@R.fun)))
        )
            warning("a covariate table has been given, yet the ",
                    sQuote("rmeasure")," function does not have ",
                    sQuote("covars")," as a formal argument: see ",
                    sQuote("?pomp"),call.=FALSE)
        if (
        (dmeasure@mode==pompfunmode$Rfun)
        &&!("covars"%in%names(formals(dmeasure@R.fun)))
        )
            warning("a covariate table has been given, yet the ",
                    sQuote("dmeasure")," function does not have ",
                    sQuote("covars")," as a formal argument: see ",
                    sQuote("?pomp"),call.=FALSE)
    }

    if ((length(tcovar)>0)&&((min(tcovar)>t0)||(max(tcovar)<max(times))))
        warning(
            "the supplied covariate covariate times ",sQuote("tcovar"),
            " do not embrace the data times: covariates may be extrapolated"
        )

    new(
        'pomp',
        rprocess = rprocess,
        dprocess = dprocess,
        dmeasure = dmeasure,
        rmeasure = rmeasure,
        dprior = dprior,
        rprior = rprior,
        skeleton = skeleton,
        skeleton.type = skeleton.type,
        skelmap.delta.t = skelmap.delta.t,
        data = data,
        times = times,
        t0 = t0,
        default.init = default.init,
        initializer = initializer,
        params = params,
        covar = covar,
        tcovar = tcovar,
        zeronames = zeronames,
        has.trans = has.trans,
        from.trans = from.trans,
        to.trans = to.trans,
        solibs = .solibs,
        userdata = userdata
    )
}

measform2pomp <- function (formulae) {
    if (!is.list(formulae))
        formulae <- list(formulae)
    nobs <- length(formulae)
    if (nobs < 1)
        stop("pomp error: to use ",sQuote("measurement.model")," you must provide at least one formula",call.=FALSE)
    for (k in seq_len(nobs)) {
        if (!inherits(formulae[[k]],"formula"))
            stop("pomp error: ",sQuote("measurement.model")," takes formulae as arguments",call.=FALSE)
    }
    obsnames <- unlist(lapply(formulae,function(x)x[[2L]]))
    distrib <- lapply(formulae,function(x)as.character(x[[3L]][[1L]]))
    ddistrib <- lapply(distrib,function(x)paste0("d",x))
    rdistrib <- lapply(distrib,function(x)paste0("r",x))
    for (k in seq_len(nobs)) {
        res <- try(
            match.fun(ddistrib[[k]]),
            silent=TRUE
        )
        if (inherits(res,'try-error'))
            stop("pomp error: distribution function ",ddistrib[[k]]," not found")
        res <- try(
            match.fun(rdistrib[[k]]),
            silent=TRUE
        )
        if (inherits(res,'try-error'))
            stop("pomp error: random deviate function ",rdistrib[[k]]," not found")
    }
    pred.args <- lapply(formulae,function(x)as.list(x[[3L]][-1L]))
    dcalls <- vector(mode='list',length=nobs)
    rcalls <- vector(mode='list',length=nobs)
    for (k in seq_len(nobs)) {
        dcalls[[k]] <- as.call(
            c(
                list(
                    as.name(ddistrib[[k]]),
                    x=obsnames[[k]]
                ),
                pred.args[[k]],
                list(
                    log=TRUE
                )
            )
        )
        rcalls[[k]] <- as.call(
            c(
                list(
                    as.name(rdistrib[[k]]),
                    n=1
                ),
                pred.args[[k]]
            )
        )
    }
    list(
        dmeasure = function (y, x, t, params, log, covars, ...) {
            f <- 0
            for (k in seq_len(nobs)) {
                f <- f+eval(
                           dcalls[[k]],
                           envir=as.list(c(y,x,params,covars,t=t))
                       )
            }
            if (log) f else exp(f)
        },
        rmeasure = function (x, t, params, covars, ...) {
            y <- numeric(length=nobs)
            names(y) <- obsnames
            for (k in seq_len(nobs)) {
                y[k] <- eval(
                    rcalls[[k]],
                    envir=as.list(c(x,params,covars,t=t))
                )
            }
            y
        }
    )
}

setMethod(
    "pomp",
    signature=signature(data="data.frame"),
    definition=function (data, times, t0, ..., rprocess, dprocess,
                         rmeasure, dmeasure, measurement.model,
                         skeleton, skeleton.type = c("map","vectorfield"),
                         skelmap.delta.t = 1,
                         initializer, rprior, dprior,
                         params, covar, tcovar,
                         obsnames, statenames, paramnames, covarnames, zeronames,
                         PACKAGE, fromEstimationScale, toEstimationScale,
                         globals, cdir, cfile) {

        data <- t(sapply(data,as.numeric))
        if ((is.numeric(times) && (times<1 || times>nrow(data))) ||
            (is.character(times) && (!(times%in%rownames(data)))) ||
            (!is.numeric(times) && !is.character(times)) ||
            length(times)!=1)
            stop("pomp error: ",sQuote("times"),
                 " must identify a single column of ",sQuote("data"))
        if (is.numeric(times)) {
            tpos <- times
        } else if (is.character(times)) {
            tpos <- match(times,rownames(data))
        }
        times <- data[tpos,,drop=TRUE]
        data <- data[-tpos,,drop=FALSE]

        pomp.constructor(
            data=data,
            times=times,
            t0=t0,
            rprocess=rprocess,
            dprocess=dprocess,
            rmeasure=rmeasure,
            dmeasure=dmeasure,
            measurement.model=measurement.model,
            dprior=dprior,
            rprior=rprior,
            skeleton=skeleton,
            skeleton.type=skeleton.type,
            skelmap.delta.t=skelmap.delta.t,
            initializer=initializer,
            params=params,
            covar=covar,
            tcovar=tcovar,
            obsnames=obsnames,
            statenames=statenames,
            paramnames=paramnames,
            covarnames=covarnames,
            zeronames=zeronames,
            PACKAGE=PACKAGE,
            fromEstimationScale=fromEstimationScale,
            toEstimationScale=toEstimationScale,
            globals=globals,
            cdir=cdir,
            cfile=cfile,
            ...
        )
    }
)

setMethod(
    "pomp",
    signature=signature(data="matrix"),
    definition=function (data, times, t0, ..., rprocess, dprocess,
                         rmeasure, dmeasure, measurement.model,
                         skeleton, skeleton.type = c("map","vectorfield"),
                         skelmap.delta.t = 1,
                         initializer, rprior, dprior, params, covar, tcovar,
                         obsnames, statenames, paramnames, covarnames, zeronames,
                         PACKAGE, fromEstimationScale, toEstimationScale,
                         globals, cdir, cfile) {

        pomp.constructor(
            data=data,
            times=times,
            t0=t0,
            rprocess=rprocess,
            dprocess=dprocess,
            rmeasure=rmeasure,
            dmeasure=dmeasure,
            measurement.model=measurement.model,
            dprior=dprior,
            rprior=rprior,
            skeleton=skeleton,
            skeleton.type=skeleton.type,
            skelmap.delta.t=skelmap.delta.t,
            initializer=initializer,
            params=params,
            covar=covar,
            tcovar=tcovar,
            obsnames=obsnames,
            statenames=statenames,
            paramnames=paramnames,
            covarnames=covarnames,
            zeronames=zeronames,
            PACKAGE=PACKAGE,
            fromEstimationScale=fromEstimationScale,
            toEstimationScale=toEstimationScale,
            globals=globals,
            cdir=cdir,
            cfile=cfile,
            ...
        )
    }
)


setMethod(
    "pomp",
    signature=signature(data="numeric"),
    definition=function (data, times, t0, ..., rprocess, dprocess,
                         rmeasure, dmeasure, measurement.model,
                         skeleton, skeleton.type = c("map","vectorfield"),
                         skelmap.delta.t = 1,
                         initializer, rprior, dprior, params, covar, tcovar,
                         obsnames, statenames, paramnames, covarnames, zeronames,
                         PACKAGE, fromEstimationScale, toEstimationScale,
                         globals, cdir, cfile) {

        pomp.constructor(
            data=matrix(data,nrow=1,ncol=length(data)),
            times=times,
            t0=t0,
            rprocess=rprocess,
            dprocess=dprocess,
            rmeasure=rmeasure,
            dmeasure=dmeasure,
            measurement.model=measurement.model,
            dprior=dprior,
            rprior=rprior,
            skeleton=skeleton,
            skeleton.type=skeleton.type,
            skelmap.delta.t=skelmap.delta.t,
            initializer=initializer,
            params=params,
            covar=covar,
            tcovar=tcovar,
            obsnames=obsnames,
            statenames=statenames,
            paramnames=paramnames,
            covarnames=covarnames,
            zeronames=zeronames,
            PACKAGE=PACKAGE,
            fromEstimationScale=fromEstimationScale,
            toEstimationScale=toEstimationScale,
            globals=globals,
            cdir=cdir,
            cfile=cfile,
            ...
        )
    }
)

setMethod(
    "pomp",
    signature=signature(data="pomp"),
    definition=function (data, times, t0, ..., rprocess, dprocess,
                         rmeasure, dmeasure, measurement.model,
                         skeleton, skeleton.type, skelmap.delta.t,
                         initializer, rprior, dprior, params, covar, tcovar,
                         obsnames, statenames, paramnames, covarnames, zeronames,
                         PACKAGE, fromEstimationScale, toEstimationScale,
                         globals, cdir, cfile) {

        if (missing(times)) times <- data@times
        if (missing(t0)) t0 <- data@t0

        mmg <- !missing(measurement.model)
        dmg <- !missing(dmeasure)
        rmg <- !missing(rmeasure)
        if (mmg) {
            if (dmg||rmg)
                warning(
                    "specifying ",sQuote("measurement.model"),
                    " overrides specification of ",
                    sQuote("rmeasure")," and ",sQuote("dmeasure")
                )
            mm <- measform2pomp(measurement.model)
            rmeasure <- mm$rmeasure
            dmeasure <- mm$dmeasure
        } else {
            if (!rmg) rmeasure <- data@rmeasure
            if (!dmg) dmeasure <- data@dmeasure
        }

        if (missing(rprocess)) rprocess <- data@rprocess
        if (missing(dprocess)) dprocess <- data@dprocess
        if (missing(rprior)) rprior <- data@rprior
        if (missing(dprior)) dprior <- data@dprior
        if (missing(initializer)) initializer <- data@initializer
        if (missing(params)) params <- data@params
        if (missing(covar)) covar <- data@covar
        if (missing(tcovar)) tcovar <- data@tcovar
        if (missing(zeronames)) zeronames <- data@zeronames
        if (missing(skeleton.type)) skeleton.type <- data@skeleton.type
        if (missing(skeleton)) skeleton <- data@skeleton
        if (missing(skelmap.delta.t)) skelmap.delta.t <- data@skelmap.delta.t
        if (missing(fromEstimationScale)) {
            if (missing(toEstimationScale)) {
                from.trans <- data@from.trans
                to.trans <- data@to.trans
            } else {
                stop("pomp error: if ",sQuote("toEstimationScale"),
                     " is supplied, then " ,
                     sQuote("fromEstimationScale")," must also be supplied")
            }
        } else {
            if (missing(toEstimationScale)) {
                stop("pomp error: if ",sQuote("fromEstimationScale"),
                     " is supplied, then " ,
                     sQuote("toEstimationScale")," must also be supplied")
            } else {
                from.trans <- fromEstimationScale
                to.trans <- toEstimationScale
            }
        }

        if (missing(obsnames)) obsnames <- character(0)
        if (missing(statenames)) statenames <- character(0)
        if (missing(paramnames)) paramnames <- character(0)
        if (missing(covarnames)) covarnames <- character(0)
        if (missing(PACKAGE)) PACKAGE <- character(0)

        pomp.constructor(
            data=data@data,
            times=times,
            t0=t0,
            rprocess=rprocess,
            dprocess=dprocess,
            rmeasure=rmeasure,
            dmeasure=dmeasure,
            dprior=dprior,
            rprior=rprior,
            skeleton=skeleton,
            skeleton.type=skeleton.type,
            skelmap.delta.t=skelmap.delta.t,
            initializer=initializer,
            covar=covar,
            tcovar=tcovar,
            obsnames=obsnames,
            statenames=statenames,
            paramnames=paramnames,
            covarnames=covarnames,
            zeronames=zeronames,
            PACKAGE=PACKAGE,
            fromEstimationScale=from.trans,
            toEstimationScale=to.trans,
            params=params,
            globals=globals,
            cdir=cdir,
            cfile=cfile,
            .solibs=data@solibs,
            userdata=data@userdata,
            ...
        )
    }
)
