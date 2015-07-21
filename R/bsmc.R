## Bayesian particle filtering codes
##
## in annotation L&W AGM == Liu & West "A General Algorithm"
## 
## params = the initial particles for the parameter values;
##          these should be drawn from the prior distribution for the parameters
## est = names of parameters to estimate; other parameters are not updated.
## smooth = parameter 'h' from AGM
## ntries = number of samplesto draw from x_{t+1} | x(k)_{t} to estimate
##          mean of mu(k)_t+1 as in sect 2.2 Liu & West
## lower  = lower bounds on prior
## upper  = upper bounds on prior

setClass(
         "bsmcd.pomp",
         contains="pomp",
         slots=c(
           transform="logical",
           post="array",
           prior="array",
           est="character",
           eff.sample.size="numeric",
           smooth="numeric",
           nfail="integer",
           cond.log.evidence="numeric",
           log.evidence="numeric"
           )
         )

bsmc.internal <- function (object, params, Np, est,
                           smooth = 0.1,
                           ntries = 1,
                           tol = 1e-17,
                           lower = -Inf, upper = Inf,
                           seed = NULL,
                           verbose = getOption("verbose"),
                           max.fail = 0,
                           transform = FALSE,
                           .getnativesymbolinfo = TRUE,
                           ...) {

  pompLoad(object)

  gnsi.rproc <- gnsi.dmeas <- as.logical(.getnativesymbolinfo)
  ptsi.inv <- ptsi.for <- TRUE
  transform <- as.logical(transform)

  if (!is.null(seed))
    warning("in ",sQuote("bsmc"),": argument ",sQuote("seed"),
            " now has no effect.  Consider using ",
            sQuote("freeze"),".")

  error.prefix <- paste(sQuote("bsmc"),"error: ")

  if (missing(params)) {
    if (length(coef(object))>0) {
      params <- coef(object)
    } else {
      stop(error.prefix,sQuote("params")," must be supplied",call.=FALSE)
    }
  }

  if (missing(Np)) Np <- NCOL(params)
  else if (is.matrix(params)&&(Np!=ncol(params)))
    warning(sQuote("Np")," is ignored when ",sQuote("params")," is a matrix")

  if ((!is.matrix(params)) && (Np > 1))
    params <- rprior(object,params=parmat(params,Np))
  
  if (transform)
    params <- partrans(object,params,dir="toEstimationScale",
                       .getnativesymbolinfo=ptsi.inv)
  ptsi.inv <- FALSE
  
  ntimes <- length(time(object))
  npars <- nrow(params)
  paramnames <- rownames(params)
  prior <- params

  if (missing(est))
    est <- paramnames[apply(params,1,function(x)diff(range(x))>0)]
  estind <- match(est,paramnames)
  npars.est <- length(estind)
  
  if (npars.est<1)
    stop(error.prefix,"no parameters to estimate",call.=FALSE)

  if (is.null(paramnames))
    stop(error.prefix,sQuote("params")," must have rownames",call.=FALSE)

  if ((length(smooth)!=1)||(smooth>1)||(smooth<=0))
    stop(error.prefix,sQuote("smooth")," must be a scalar in [0,1)",call.=FALSE)

  hsq <- smooth^2             #  see Liu & West eq(3.6) p10
  shrink <- sqrt(1-hsq)

  if (
      ((length(lower)>1)&&(length(lower)!=npars.est))||
      ((length(upper)>1)&&(length(upper)!=npars.est))
      ) {
    stop(
         error.prefix,
         sQuote("lower")," and ",sQuote("upper"),
         " must each have length 1 or length equal to that of ",sQuote("est"),
         call.=FALSE                   
         )
  }

  for (j in seq_len(Np)) {
    if (any((params[estind,j]<lower)|(params[estind,j]>upper))) {
      ind <- which((params[estind,j]<lower)|(params[estind,j]>upper))
      stop(
           error.prefix,
           "parameter(s) ",paste(paramnames[estind[ind]],collapse=","),
           " in column ",j," in ",sQuote("params"),
           " is/are outside the box defined by ",
           sQuote("lower")," and ",sQuote("upper"),
           call.=FALSE
           )
    }
  }

  xstart <- init.state(
                       object,
                       params=if (transform) {
                         partrans(object,params,dir="fromEstimationScale",
                                  .getnativesymbolinfo=ptsi.for)
                       } else {
                         params
                       }
                       )
  statenames <- rownames(xstart)
  nvars <- nrow(xstart)
  ptsi.for <- FALSE
  
  times <- time(object,t0=TRUE)
  x <- xstart

  evidence <- rep(NA,ntimes)
  eff.sample.size <- rep(NA,ntimes)
  nfail <- 0
  
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
      cat("at step",nt,"(time =",times[nt+1],")\n")
      print(
            rbind(
                  prior.mean=params.mean[estind],
                  prior.sd=sqrt(diag(params.var))
                  )
            )
    }

    ## update mean of states at time nt as per L&W AGM (1) 
    tries <- rprocess(
                      object,
                      xstart=parmat(x,nrep=ntries),
                      times=times[c(nt,nt+1)],
                      params=if (transform) {
                        partrans(object,params,dir="fromEstimationScale",
                                 .getnativesymbolinfo=ptsi.for)
                      } else {
                        params
                      },
                      offset=1,
                      .getnativesymbolinfo=gnsi.rproc
                      )
    dim(tries) <- c(nvars,Np,ntries,1)
    mu <- apply(tries,c(1,2,4),mean)
    rownames(mu) <- statenames
    ## shrink parameters towards mean as per Liu & West eq (3.3) and L&W AGM (1)
    m <- shrink*params+(1-shrink)*params.mean
    gnsi.rproc <- FALSE
    
    ## evaluate probability of obervation given mean value of parameters and states (used in L&W AGM (5) below)
    g <- dmeasure( 
                  object,
                  y=object@data[,nt,drop=FALSE],
                  x=mu,
                  times=times[nt+1],
                  params=if (transform) {
                    partrans(object,m,dir="fromEstimationScale",
                             .getnativesymbolinfo=ptsi.for)
                  } else {
                    m
                  },
                  .getnativesymbolinfo=gnsi.dmeas
                  )	
    gnsi.dmeas <- FALSE
    storeForEvidence1 <- log(sum(g))
    ## sample indices -- From L&W AGM (2)
    ##              k <- .Call(systematic_resampling,g)
    k <- sample.int(n=Np,size=Np,replace=TRUE,prob=g)
    params <- params[,k]
    m <- m[,k]
    g <- g[k]

    ## sample new parameter vector as per L&W AGM (3) and Liu & West eq(3.2)
    pvec <- try(
                rmvnorm(
                        n=Np,
                        mean=rep(0,npars.est),
                        sigma=hsq*params.var,
                        method="svd"
                        ),
                silent=FALSE
                )
    if (inherits(pvec,"try-error"))
      stop(error.prefix,"error in ",sQuote("rmvnorm"),call.=FALSE)
    if (!all(is.finite(pvec)))
      stop(error.prefix,"extreme particle depletion",call.=FALSE)
    params[estind,] <- m[estind,]+t(pvec)

    if (transform)
      tparams <- partrans(object,params,dir="fromEstimationScale",
                          .getnativesymbolinfo=ptsi.for)
    
    ## sample current state vector x^(g)_(t+1) as per L&W AGM (4)
    X <- rprocess(
                  object,
                  xstart=x[,k,drop=FALSE],
                  times=times[c(nt,nt+1)],
                  params=if (transform) {
                    tparams
                  } else {
                    params
                  },
                  offset=1,
                  .getnativesymbolinfo=gnsi.rproc
                  )

    ## evaluate likelihood of observation given X (from L&W AGM (4))
    numer <- dmeasure(
                      object,
                      y=object@data[,nt,drop=FALSE],
                      x=X,
                      times=times[nt+1],
                      params=if (transform) {
                        tparams
                      } else {
                        params
                      },
                      .getnativesymbolinfo=gnsi.dmeas
                      )
    ## evaluate weights as per L&W AGM (5)

    weights <- numer/g
    storeForEvidence2 <- log(mean(weights))
    
    ## apply box constraints as per the priors          
    for (j in seq_len(Np)) {
      ## the following seems problematic: will it tend to make the boundaries repellors
      if (any((params[estind,j]<lower)|(params[estind,j]>upper))) {
        weights[j] <- 0 
      }
      ## might this rejection method be preferable?
      ## while (any((params[estind,j]<lower)|(params[estind,j]>upper))) {
      ##   ## rejection method
      ##   pvec <- try(
      ##               mvtnorm::rmvnorm(
      ##                                n=1,
      ##                                mean=rep(0,npars.est),
      ##                                sigma=hsq*params.var,
      ##                                method="eigen"
      ##                                ),
      ##               silent=FALSE
      ##               )
      ##   if (inherits(pvec,"try-error"))
      ##     stop(error.prefix,"error in ",sQuote("rmvnorm"),call.=FALSE)
      ##   if (!all(is.finite(pvec)))
      ##     stop(error.prefix,"extreme particle depletion",call.=FALSE)
      ##   params[estind,j] <- m[estind,j]+pvec[1,]
      ## }
    }

    x[,] <- X                
    
    ## test for failure to filter
    dim(weights) <- NULL
    failures <- ((weights<tol)|(!is.finite(weights))) # test for NA weights
    all.fail <- all(failures)
    if (all.fail) {                     # all particles are lost
      if (verbose) {
        message("filtering failure at time t = ",times[nt+1])
      }
      nfail <- nfail+1
      if (nfail > max.fail)
        stop(error.prefix,"too many filtering failures",call.=FALSE)
      evidence[nt] <- log(tol)          # worst log-likelihood
      weights <- rep(1/Np,Np)
      eff.sample.size[nt] <- 0
    } else {                  # not all particles are lost
      ## compute log-likelihood
      evidence[nt] <- storeForEvidence1+storeForEvidence2
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
      ## smp <- .Call(systematic_resampling,weights)
      smp <- sample.int(n=Np,size=Np,replace=TRUE,prob=weights)
      x <- x[,smp,drop=FALSE]
      params[estind,] <- params[estind,smp,drop=FALSE]
    }

    .getnativesymbolinfo <- FALSE
    
  }

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
          "bsmc",
          signature=signature(object="pomp"),
          definition = function (object, params, Np, est,
            smooth = 0.1,
            ntries = 1,
            tol = 1e-17,
            lower = -Inf, upper = Inf,
            verbose = getOption("verbose"),
            max.fail = 0,
            transform = FALSE,
            ...) {
            bsmc.internal(
                          object=object,
                          params=params,
                          Np=Np,
                          est=est,
                          smooth=smooth,
                          ntries=ntries,
                          tol=tol,
                          lower=lower,
                          upper=upper,
                          verbose=verbose,
                          max.fail=max.fail,
                          transform=transform,
                          ...
                          )
          }
          )

setMethod("$",
          signature(x="bsmcd.pomp"),
          definition = function (x, name) slot(x,name)
          )

bsmc.plot <- function (prior, post, pars, thin, ...) {
  p1 <- sample.int(n=ncol(prior),size=min(thin,ncol(prior)))
  p2 <- sample.int(n=ncol(post),size=min(thin,ncol(post)))
  if (!all(pars%in%rownames(prior))) {
    missing <- which(!(pars%in%rownames(prior)))
    stop("unrecognized parameters: ",paste(sQuote(pars[missing]),collapse=","))
    
  }
  prior <- t(prior[pars,])
  post <- t(post[pars,])
  all <- rbind(prior,post)
  pairs(
        all,
        labels=pars,
        panel=function (x, y, ...) { ## prior, posterior pairwise scatterplot
          op <- par(new=TRUE)
          on.exit(par(op))
          i <- which(x[1L]==all[1L,])
          j <- which(y[1L]==all[1L,])
          points(prior[p1,i],prior[p1,j],pch=20,col=rgb(0.85,0.85,0.85,0.1),xlim=range(all[,i]),ylim=range(all[,j]))
          points(post[p2,i],post[p2,j],pch=20,col=rgb(0,0,1,0.01))
        },
        diag.panel=function (x, ...) { ## marginal posterior histogram
          i <- which(x[1L]==all[1L,])
          d1 <- density(prior[,i])
          d2 <- density(post[,i])
          usr <- par('usr')
          op <- par(usr=c(usr[c(1L,2L)],0,1.5*max(d1$y,d2$y)))
          on.exit(par(op))
          polygon(d1,col=rgb(0.85,0.85,0.85,0.5))
          polygon(d2,col=rgb(0,0,1,0.5))
        }
        )
}

setMethod(
          "plot",
          signature(x="bsmcd.pomp"),
          function (x, y, ..., pars, thin) {
            if (!missing(y))
              warning(sQuote("y")," is ignored")
            if (missing(pars)) pars <- x@est
            if (missing(thin)) thin <- Inf
            bsmc.plot(
                      prior=if (x@transform) partrans(x,x@prior,dir="fromEstimationScale") else x@prior,
                      post=if (x@transform) partrans(x,x@post,dir="fromEstimationScale") else x@post,
                      pars=pars,
                      thin=thin,
                      ...
                      )
          }
          )
