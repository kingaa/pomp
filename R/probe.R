setClass(
    "probed.pomp",
    contains="pomp",
    slots=c(
        probes="list",
        datvals="numeric",
        simvals="array",
        quantiles="numeric",
        pvals="numeric",
        synth.loglik="numeric",
        seed="integer"
    )
)

probe.internal <- function (object, probes, params, nsim = 1, seed = NULL, ...) {

    pompLoad(object)

    if (!is.list(probes)) probes <- list(probes)
    if (!all(sapply(probes,is.function)))
        stop(sQuote("probes")," must be a function or a list of functions")
    if (!all(sapply(probes,function(f)length(formals(f))==1)))
        stop("each probe must be a function of a single argument")

    seed <- as.integer(seed)
    
    if (missing(params)) params <- coef(object)

    ## apply probes to data
    datval <- .Call(apply_probe_data,object,probes)
    nprobes <- length(datval)
    if (nprobes > nsim)
        stop(sQuote("nsim"),"(=",nsim,") should be (much) larger than the number of probes (=",nprobes,")")

    ## apply probes to model simulations
    simval <- .Call(
        apply_probe_sim,
        object=object,
        nsim=nsim,
        params=params,
        seed=seed,
        probes=probes,
        datval=datval
    )
    
    pvals <- numeric(nprobes)
    names(pvals) <- names(datval)
    quants <- numeric(nprobes)
    names(quants) <- names(datval)
    for (k in seq_len(nprobes)) {
        r <- min(sum(simval[,k]>datval[k]),sum(simval[,k]<datval[k]))
        tails <- (r+1)/(nsim+1)
        pvals[k] <- min(2*tails,1)
        quants[k] <- sum(simval[,k]<datval[k])/nsim
    }

    ll <- .Call(synth_loglik,simval,datval)

    coef(object) <- params

    pompUnload(object)

    new(
        "probed.pomp",
        object,
        probes=probes,
        datvals=datval,
        simvals=simval,
        quantiles=quants,
        pvals=pvals,
        synth.loglik=ll,
        seed=seed
    )
}

setMethod(
    "probe",
    signature=signature(object="pomp"),
    definition=function (object, probes, params, nsim = 1, seed = NULL, ...)
    {
        probe.internal(
            object=object,
            probes=probes,
            params=params,
            nsim=nsim,
            seed=seed,
            ...
        )
    }
)

setMethod(
    "probe",
    signature=signature(object="probed.pomp"),
    definition=function (object, probes, params, nsim, seed, ...) {

        if (missing(probes)) probes <- object@probes
        if (missing(nsim)) nsim <- nrow(object@simvals)
        if (missing(seed)) seed <- object@seed

        probe(
            object=as(object,"pomp"),
            probes=probes,
            nsim=nsim,
            seed=seed,
            ...
        )
    }
)

probeplot.internal <- function (x, ...) {
    ##function for plotting diagonal panels
    diag.panel.hist <- function(x, ...) {
        ##plot a histogram for the simulations
        usr <- par("usr")
        on.exit(par(usr))
        par(usr=c(usr[c(1L,2L)],0,1.5))
        h <- hist(x[-1L],plot=FALSE)
        breaks <- h$breaks
        nB <- length(breaks)
        y <- h$counts
        y <- y/max(y)
        rect(breaks[-nB],0,breaks[-1L],y,...)
        ##plot the data point
        lines(c(x[1L],x[1L]),c(0,max(h$counts)),col="red")
    }

    ##function for plotting above-diagonal panels
    above.diag.panel <- function (x, y, ...) {
        ##plot the simulations
        points(x[-1L],y[-1L],...)
        ##plot the data
        mMx <- c(min(x),max(x))
        mMy <- c(min(y),max(y))
        lines(c(x[1L],x[1L]),mMy,col="red")
        lines(mMx,c(y[1L],y[1L]),col="red")
    }
    
    ##function for plotting below-diagonal panels
    below.diag.panel <- function (x, y, ...) {
        mMx <- c(min(x),max(x))
        mMy <- c(min(y),max(y))
        x <- x[-1L]
        y <- y[-1L]
        correls <- round(cor(x,y),3)
        text(mean(mMx),mean(mMy),correls,cex=1)
    }
    
    ##prepare the arguments for these functions
    nprobes <- length(x@datvals)
    nsim <- nrow(x@simvals)
    datsimvals <- array(dim=c(nsim+1,nprobes))
    datsimvals[1L,] <- x@datvals
    datsimvals[-1L,] <- x@simvals
    
    labels <- paste("pb",seq_len(nprobes))
    if (!is.null(names(x@datvals)))
        labels <- ifelse(names(x@datvals)=="",labels,names(x@datvals))
    lab.plus <- paste(labels,paste0("p=",round(x@pvals,3)),sep="\n")
    ##now make the plot

    if (nprobes>1) {
        pairs(
            datsimvals,
            diag.panel=diag.panel.hist,
            lower.panel=below.diag.panel,
            upper.panel=above.diag.panel,
            labels=lab.plus,
            cex.labels=if (nprobes>5) 5/nprobes else 1
        )
    } else {
        plot(datsimvals,datsimvals,type="n",xlab="",ylab="",yaxt="n",main=lab.plus)
        diag.panel.hist(datsimvals)
    }
}

setMethod("plot",
          signature=signature(x="probed.pomp"), 
          definition=function (x, y, ...) {
              if (!missing(y))
                  warning(sQuote("y")," is ignored")
              probeplot.internal(x=x,...)
          }
          )

setMethod(
    "summary",
    signature=signature(object="probed.pomp"),
    definition=function (object, ...) {
        list(
            coef=coef(object),
            nsim=nrow(object@simvals),
            quantiles=object@quantiles,
            pvals=object@pvals,
            synth.loglik=object@synth.loglik
        )
    }
)

setAs(
    from="probed.pomp",
    to="data.frame",
    def = function (from) {
        x <- rbind(from@datvals,as.data.frame(from@simvals))
        rownames(x) <- c(
            "data",
            paste("sim",seq_len(nrow(from@simvals)),sep=".")
        )
        x
    }
)

as.data.frame.probed.pomp <- function (x, row.names, optional, ...) as(x,"data.frame")

setMethod(
    "logLik",
    signature=signature(object="probed.pomp"),
    definition=function(object,...)object@synth.loglik
)

setMethod(
    "$",
    signature=signature(x="probed.pomp"),
    definition=function(x, name)slot(x,name)
)

values.probe.internal <- function (object, ...) {
    x <- as.data.frame(rbind(object@datvals,object@simvals))
    row.names(x) <- seq.int(from=0,to=nrow(x)-1)
    x$.id <- factor(c("data",rep("sim",nrow(x)-1)))
    x
}

setMethod(
    "values",
    signature=signature(object="probed.pomp"),
    definition=values.probe.internal
)
