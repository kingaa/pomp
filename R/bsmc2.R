## Bayesian particle filtering codes
##
## in annotation L&W AGM == Liu & West "A General Algorithm"
## 
## params = the initial particles for the parameter values;
##          these should be drawn from the prior distribution for the parameters
## est = names of parameters to estimate; other parameters are not updated.
## smooth = parameter 'h' from AGM

bsmc2.internal <- function (object, params, Np, est,
                            smooth, tol, seed = NULL,
                            verbose = getOption("verbose"),
                            max.fail, transform, .getnativesymbolinfo = TRUE,
                            ...) {

    ep <- paste0("in ",sQuote("bsmc2"),": ")

    object <- as(object,"pomp")
    pompLoad(object)
    
    gnsi.rproc <- gnsi.dmeas <- as.logical(.getnativesymbolinfo)
    ptsi.inv <- ptsi.for <- TRUE
    transform <- as.logical(transform)

    if (!is.null(seed))
        warning(ep,"argument ",sQuote("seed"),
                " now has no effect.  Consider using ",
                sQuote("freeze"),".",call.=FALSE)

    if (missing(params)) {
        if (length(coef(object))>0) {
            params <- coef(object)
        } else {
            stop(ep,sQuote("params")," must be supplied",call.=FALSE)
        }
    }

    if (missing(Np)) Np <- NCOL(params)
    else if (is.matrix(params)&&(Np!=ncol(params)))
        warning(ep,sQuote("Np")," is ignored when ",sQuote("params"),
                " is a matrix",call.=FALSE)

    if ((!is.matrix(params)) && (Np > 1)) {
        params <- tryCatch(
            rprior(object,params=parmat(params,Np)),
            error = function (e) {
                stop(ep,sQuote("rprior")," error: ",
                     conditionMessage(e),call.=FALSE)
            }
        )
    }
    
    if (transform)
        params <- partrans(object,params,dir="toEstimationScale",
                           .getnativesymbolinfo=ptsi.inv)
    ptsi.inv <- FALSE
    
    params <- as.matrix(params)

    ntimes <- length(time(object))
    npars <- nrow(params)
    paramnames <- rownames(params)
    prior <- params

    if (missing(est))
        est <- paramnames[apply(params,1,function(x)diff(range(x))>0)]
    estind <- match(est,paramnames)
    npars.est <- length(estind)
    
    if (npars.est<1)
        stop(ep,"no parameters to estimate",call.=FALSE)

    if (is.null(paramnames))
        stop(ep,sQuote("params")," must have rownames",call.=FALSE)

    if ((length(smooth)!=1)||(smooth>1)||(smooth<=0))
        stop(ep,sQuote("smooth")," must be a scalar in [0,1)",call.=FALSE)

    hsq <- smooth^2             #  see Liu & West eq(10.3.12)
    shrink <- sqrt(1-hsq)       #  'a' parameter of Liu & West

    xstart <- init.state(
        object,
        params=if (transform) {
                   partrans(object,params,dir="fromEstimationScale",
                            .getnativesymbolinfo=ptsi.for)
               } else {
                   params
               }
    )
    nvars <- nrow(xstart)
    ptsi.for <- FALSE
    
    times <- time(object,t0=TRUE)
    x <- xstart

    evidence <- as.numeric(rep(NA,ntimes))
    eff.sample.size <- as.numeric(rep(NA,ntimes))
    nfail <- 0L
    
    mu <- array(data=NA,dim=c(nvars,Np,1)) 
    rownames(mu) <- rownames(xstart)
    m  <- array(data=NA,dim=c(npars,Np))
    rownames(m) <- rownames(params)
    
    for (nt in seq_len(ntimes)) {
        
        ## calculate particle means ; as per L&W AGM (1)
        params.mean <- apply(params,1,mean) 
        ## calculate particle covariances : as per L&W AGM (1)
        params.var  <- cov(t(params[estind,,drop=FALSE]))

        if (verbose) {
            cat("at step ",nt," (time =",times[nt+1],")\n",sep="")
            print(
                rbind(
                    prior.mean=params.mean[estind],
                    prior.sd=sqrt(diag(params.var))
                )
            )
        }

        m <- shrink*params+(1-shrink)*params.mean
        
        ## sample new parameter vector as per L&W AGM (3) and Liu & West eq(3.2)
        pert <- tryCatch(
            rmvnorm(
                n=Np,
                mean=rep(0,npars.est),
                sigma=hsq*params.var,
                method="svd"
            ),
            error = function (e) {
                stop(ep,"in ",sQuote("rmvnorm"),": ",conditionMessage(e),call.=FALSE)
            }
        )
        if (!all(is.finite(pert)))
            stop(ep,"extreme particle depletion",call.=FALSE)

        params[estind,] <- m[estind,]+t(pert)

        if (transform)
            tparams <- partrans(object,params,dir="fromEstimationScale",
                                .getnativesymbolinfo=ptsi.for)
        
        xpred <- rprocess(
            object,
            xstart=x,
            times=times[c(nt,nt+1)],
            params=if (transform) {
                       tparams
                   } else {
                       params
                   },
            offset=1,
            .getnativesymbolinfo=gnsi.rproc
        )
        gnsi.rproc <- FALSE

        ## evaluate likelihood of observation given xpred (from L&W AGM (4))
        weights <- tryCatch(
            dmeasure(
                object,
                y=object@data[,nt,drop=FALSE],
                x=xpred,
                times=times[nt+1],
                params=if (transform) {
                           tparams
                       } else {
                           params
                       },
                .getnativesymbolinfo=gnsi.dmeas
            ),
            error = function (e) {
                stop(ep,sQuote("dmeasure")," error: ",conditionMessage(e),
                     call.=FALSE)
            }
        )
            
        gnsi.dmeas <- FALSE
        ## evaluate weights as per L&W AGM (5)

        storeForEvidence <- log(mean(weights))
        
        x[,] <- xpred
        
        ## test for failure to filter
        dim(weights) <- NULL   ### needed?
        failures <- ((weights<tol)|(!is.finite(weights))) # test for NA weights
        all.fail <- all(failures)
        if (all.fail) {                     # all particles are lost
            if (verbose) {
                message("filtering failure at time t = ",times[nt+1])
            }
            nfail <- nfail+1
            if (nfail > max.fail)
                stop(ep,"too many filtering failures",call.=FALSE)
            evidence[nt] <- log(tol)          # worst log-likelihood
            weights <- rep(1/Np,Np)
            eff.sample.size[nt] <- 0
        } else {                  # not all particles are lost
            ## compute log-likelihood
            evidence[nt] <- storeForEvidence
            weights[failures] <- 0
            weights <- weights/sum(weights)
            ## compute effective sample-size
            eff.sample.size[nt] <- 1/crossprod(weights)
        }

        if (verbose) {
            cat("effective sample size =",round(eff.sample.size[nt],1),"\n")
        }

        ## Matrix with samples (columns) from filtering distribution theta.t | Y.t
        if (!all.fail) {
            smp <- tryCatch(
                .Call(systematic_resampling,weights),
                error = function (e) {
                    stop(ep,conditionMessage(e),call.=FALSE)
                }
            )
            x <- x[,smp,drop=FALSE]
            params[estind,] <- params[estind,smp,drop=FALSE]
        }

        .getnativesymbolinfo <- FALSE
        
    }
    
    if (nfail>0)
        warning(
            ep,nfail,
            ngettext(
                nfail,
                msg1=" filtering failure occurred.",
                msg2=" filtering failures occurred."
            ),
            call.=FALSE
        )

    ## replace parameters with point estimate (posterior median)
    coef(object,transform=transform) <- apply(params,1,median)

    pompUnload(object)

    new(
        "bsmcd.pomp",
        object,
        transform=transform,
        post=params,
        prior=prior,
        est=as.character(est),
        eff.sample.size=eff.sample.size,
        smooth=smooth,
        nfail=as.integer(nfail),
        cond.log.evidence=evidence,
        log.evidence=sum(evidence)
    )
}

setMethod(
    "bsmc2",
    signature=signature(object="pomp"),
    definition = function (object, params, Np, est,
                           smooth = 0.1, tol = 1e-17,
                           verbose = getOption("verbose"),
                           max.fail = 0, transform = FALSE,
                           ...) {
        bsmc2.internal(
            object=object,
            params=params,
            Np=Np,
            est=est,
            smooth=smooth,
            tol=tol,
            verbose=verbose,
            max.fail=max.fail,
            transform=transform,
            ...
        )
    }
)
