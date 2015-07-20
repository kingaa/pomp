# Authors:
# Cai GoGwilt, MIT
# Daniel Reuman, Imperial College London
# Aaron A. King, U of Michigan

setClass(
         "spect.pomp",
         contains="pomp",
         slots=c(
           kernel.width="numeric",
           transform="function",
           freq="numeric",
           datspec="array",
           simspec="array",
           pvals="numeric",
           detrend="character"
           )
         )

setMethod(
          "summary",
          "spect.pomp",
          function (object, ...) {
            list(
                 coef=coef(object),
                 nsim=nrow(object@simspec),
                 pvals=object@pvals
                 )
          }
          )

## detrends in one of several ways, according to type.
## tseries is a numeric vector,
pomp.detrend <- function (tseries, type)
  switch(
         type,
         mean=tseries-mean(tseries),
         linear={
           m <- cbind(1,seq_along(tseries))
           .lm.fit(m,tseries)$residuals
         },
         quadratic={
           x <- seq_along(tseries)
           m <- cbind(1,x,x*x)
           .lm.fit(m,tseries)$residuals
         },
         tseries
         )

## The default smoothing kernel for the R spec.pgram function is weird.
## This function creates a better one.
reuman.kernel <- function (kernel.width) {
  ker <- kernel("modified.daniell",m=kernel.width)
  x <- seq.int(from=0,to=kernel.width,by=1)/kernel.width
  ker[[1L]] <- (15/(16*2*pi))*((x-1)^2)*((x+1)^2)
  ker[[1L]] <- ker[[1L]]/(2*sum(ker[[1L]][-1])+ker[[1L]][1L])
  attr(ker,"name") <- NULL
  ker
}

compute.spect.data <- function (object, vars, transform, detrend, ker) {
  dat <- obs(object,vars)
  if (any(is.na(dat)))
    stop(sQuote("spect")," is incompatible with NAs in the data")
  dt <- diff(time(object,t0=FALSE))
  base.freq <- 1/mean(dt)
  dt.tol <- 0.025
  if (max(dt)-min(dt)>dt.tol*mean(dt))
    stop(sQuote("spect")," assumes evenly spaced times")
  for (j in seq_along(vars)) {
    sp <- spec.pgram(
                     pomp.detrend(
                                  transform(dat[j,]),
                                  type=detrend
                                  ),
                     spans=ker,
                     taper=0,
                     pad=0,
                     fast=FALSE,
                     detrend=FALSE,
                     plot=FALSE
                     )
    if (j==1) {
      freq <- base.freq*sp$freq
      datspec <- array(
                       dim=c(length(freq),nrow(dat)),
                       dimnames=list(NULL,vars)
                       )
    }
    datspec[,j] <- log10(sp$spec)
  }
  list(freq=freq,spec=datspec)
}

compute.spect.sim <- function (object, params, vars, nsim, seed, transform, detrend, ker) {
  sims <- try(
              simulate(
                       object,
                       nsim=nsim,
                       seed=seed,
                       params=params,
                       obs=TRUE,
                       times=time(object,t0=FALSE),
                       t0=timezero(object)
                       ),
              silent=FALSE
              )
  if (inherits(sims,"try-error"))
    stop(sQuote("spect")," error: cannot simulate")
  sims <- sims[vars,,,drop=FALSE]
  if (any(is.na(sims)))
    stop("NA in simulated data series")
  nobs <- length(vars)
  for (j in seq_len(nobs)) {
    for (k in seq_len(nsim)) {
      sp <- spec.pgram(
                       pomp.detrend(
                                    transform(sims[j,k,]),
                                    type=detrend
                                    ),
                       spans=ker,
                       taper=0,
                       pad=0,
                       fast=FALSE,
                       detrend=FALSE,
                       plot=FALSE
                       )
      if ((j==1)&&(k==1)) {
        simspec <- array(
                         dim=c(nsim,length(sp$freq),nobs),
                         dimnames=list(NULL,NULL,vars)
                         )
      }
      simspec[k,,j] <- log10(sp$spec)
    }
  }
  simspec
}

setMethod(
          "spect",
          signature(object="pomp"),
          function (object, params, vars, kernel.width, nsim, seed = NULL,
                    transform = identity,
                    detrend = c("none","mean","linear","quadratic"),
                    ...) {

            pompLoad(object)
            
            if (missing(params)) params <- coef(object)

            if (missing(vars))
              vars <- rownames(object@data)
            
            if (missing(kernel.width))
              stop(sQuote("kernel.width")," must be specified")
            if (missing(nsim)||(nsim<1))
              stop(sQuote("nsim")," must be specified as a positive integer")

            detrend <- match.arg(detrend)
            ker <- reuman.kernel(kernel.width)

            ds <- compute.spect.data(
                                     object,
                                     vars=vars,
                                     transform=transform,
                                     detrend=detrend,
                                     ker=ker
                                     )
            freq <- ds$freq
            datspec <- ds$spec
            simspec <- compute.spect.sim(
                                         object,
                                         params=params,
                                         vars=vars,
                                         nsim=nsim,
                                         seed=seed,
                                         transform=transform,
                                         detrend=detrend,
                                         ker=ker
                                         )

            pvals <- numeric(length(vars)+1)
            names(pvals) <- c(vars,"all")
            mean.simspec <- colMeans(simspec) # mean spectrum of simulations
            totdatdist <- 0
            totsimdist <- 0
            for (j in seq_along(vars)) {
              ## L-2 distance between data and mean simulated spectrum
              datdist <- sum((datspec[,j]-mean.simspec[,j])^2)
              ## L-2 distance betw. each sim. and mean simulated spectrum
              simdist <- vapply(
                                seq_len(nsim),
                                function(k)sum((simspec[k,,j]-mean.simspec[,j])^2),
                                numeric(1)
                                )
              pvals[j] <- (nsim+1-sum(simdist<datdist))/(nsim+1)
              totdatdist <- totdatdist+datdist
              totsimdist <- totsimdist+simdist
            }
            pvals[length(vars)+1] <- (nsim+1-sum(totsimdist<totdatdist))/(nsim+1)

            coef(object) <- params

            pompUnload(object)
            
            new(
                "spect.pomp",
                object,
                kernel.width=kernel.width,
                transform=transform,
                detrend=detrend,
                freq=freq,
                datspec=datspec,
                simspec=simspec,
                pvals=pvals
                )
          }
          )

setMethod(
          "spect",
          signature(object="spect.pomp"),
          function (object, params, vars, kernel.width, nsim, seed = NULL, transform, detrend, ...) {
            if (missing(params)) params <- coef(object)
            if (missing(vars)) vars <- rownames(object@datspec)
            if (missing(kernel.width)) kernel.width <- object@kernel.width
            if (missing(nsim)) nsim <- nrow(object@simspec)
            if (missing(transform)) transform <- object@transform
            if (missing(detrend)) detrend <- object@detrend
            spect(
                  as(object,"pomp"),
                  params=params,
                  vars=vars,
                  kernel.width=kernel.width,
                  nsim=nsim,
                  seed=seed,
                  transform=transform,
                  detrend=detrend,
                  ...
                  )
          }
          )

setMethod(
          "plot",
          "spect.pomp",
          function (x, y, max.plots.per.page = 4,
                    plot.data = TRUE,
                    quantiles = c(.025, .25, .5, .75, .975),
                    quantile.styles = list(lwd=1, lty=1, col="gray70"),
                    data.styles = list(lwd=2, lty=2, col="black")) {
            spomp <- x
            nquants <- length(quantiles)
            
            if (is.list(quantile.styles)) {
              for (i in c("lwd", "lty", "col")) {
                if (is.null(quantile.styles[[i]]))
                  quantile.styles[[i]] <- rep(1,nquants)
                if (length(quantile.styles[[i]])==1)
                  quantile.styles[[i]] <- rep(quantile.styles[[i]],nquants)
                if (length(quantile.styles[[i]])<nquants) {
                  warning(
                          sQuote("quantile.styles"),
                          " contains an element with more than 1 entry but fewer entries than quantiles"
                          )
                  quantile.styles[[i]]<-rep(quantile.styles[[i]],nquants)
                }
              }
            } else {
              stop(sQuote("quantile.styles")," must be a list")
            }
            
            if (plot.data) {
              nreps <- ncol(spomp@datspec)
              if (is.list(data.styles)) {
                for (i in c("lwd", "lty", "col")) {
                  if(is.null(data.styles[[i]]))
                    data.styles[[i]] <- rep(2,nreps)
                  if(length(data.styles[[i]])==1)
                    data.styles[[i]] <- rep(data.styles[[i]],nreps)
                  if(length(data.styles[[i]]) < nreps) {
                    warning(
                            sQuote("data.styles"),
                            "contains an element with more than 1 entry but fewer entries of observed variables"
                            )
                    data.styles[[i]] <- rep(data.styles[[i]],nreps)
                  }
                }
              } else {
                stop(sQuote("data.styles")," must be a list")
              }
            }

            dimsim <- dim(spomp@simspec)
            nfreq <- dimsim[2]
            nobs <- dimsim[3]
            oldpar <- par(
                          mfrow=c(min(nobs,max.plots.per.page),1),
                          mar=c(3,3,1,0.5),
                          mgp=c(2,1,0),
                          bty="l",
                          ask=if (nobs>max.plots.per.page) TRUE else par("ask")
                          )
            on.exit(par(oldpar))
            ylabs <- dimnames(spomp@simspec)[[3]]
            for (i in seq_len(nobs)) {
              spectraquants <- array(dim=c(nfreq,length(quantiles)))
              for (j in seq_len(nfreq))
                spectraquants[j,] <- quantile(
                                              spomp@simspec[,j,i],
                                              probs=quantiles
                                              )
              if (plot.data) {
                ylimits <- c(
                             min(spectraquants,spomp@datspec[,i]),
                             max(spectraquants,spomp@datspec[,i])
                             )
              } else {
                ylimits <- c(
                             min(spectraquants),
                             max(spectraquants)
                             )
              }
              plot(
                   NULL,
                   xlim=range(spomp@freq),ylim=ylimits,
                   xlab=if (i==nobs) "frequency" else "",
                   ylab=expression(paste(log[10],"power"))
                   )
              title(
                    main=paste0(ylabs[i],", p = ",round(spomp@pvals[i],4)),
                    line=0
                    )
              for (j in seq_along(quantiles)) {
                lines(
                      x=spomp@freq,
                      y=spectraquants[,j],
                      lwd=quantile.styles$lwd[j], 
                      lty=quantile.styles$lty[j],
                      col=quantile.styles$col[j]
                      )
              }
              if(plot.data) {
                lines(
                      x=spomp@freq,
                      y=spomp@datspec[,i],
                      lty=data.styles$lty[i], 
                      lwd=data.styles$lwd[i],
                      col=data.styles$col[i]
                      )
              }
            }
          }
          )

