## this file defines methods for the 'pmcmc' and 'pmcmcList' classes

## extract the estimated log likelihood
setMethod('logLik','pmcmc',function(object,...)object@loglik)

## pmcmcList class
setClass(
         'pmcmcList',
         contains='list',
         validity=function (object) {
           if (!all(sapply(object,is,'pmcmc'))) {
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
                              ": to be combined, ",sQuote("pmcmc"),
                              " objects must have chains of equal length"
                              )
             return(retval)
           }
           TRUE
         }
         )

setMethod(
          'c',
          signature=signature(x='pmcmc'),
          definition=function (x, ...) {
            y <- list(...)
            if (length(y)==0) {
              new("pmcmcList",list(x))
            } else {
              p <- sapply(y,is,'pmcmc')
              pl <- sapply(y,is,'pmcmcList')
              if (!all(p||pl))
                stop("cannot mix ",sQuote("pmcmc"),
                     " and non-",sQuote("pmcmc")," objects")
              y[p] <- lapply(y[p],list)
              y[pl] <- lapply(y[pl],as,"list")
              new("pmcmcList",c(list(x),y,recursive=TRUE))
            }
          }
          )

setMethod(
          'c',
          signature=signature(x='pmcmcList'),
          definition=function (x, ...) {
            y <- list(...)
            if (length(y)==0) {
              x
            } else {
              p <- sapply(y,is,'pmcmc')
              pl <- sapply(y,is,'pmcmcList')
              if (!all(p||pl))
                stop("cannot mix ",sQuote("pmcmc"),
                     " and non-",sQuote("pmcmc")," objects")
              y[p] <- lapply(y[p],list)
              y[pl] <- lapply(y[pl],as,"list")
              new("pmcmcList",c(as(x,"list"),y,recursive=TRUE))
            }
          }
          )

setMethod(
          "[",
          signature=signature(x="pmcmcList"),
          definition=function(x, i, ...) {
            new('pmcmcList',as(x,"list")[i])
          }
          )

## extract the convergence record as a coda::mcmc object
setMethod(
          'conv.rec',
          signature=signature(object='pmcmc'),
          function (object, pars, ...) {
            if (missing(pars)) pars <- colnames(object@conv.rec)
            coda::mcmc(object@conv.rec[,pars,drop=FALSE])
          }
          )

## extract the convergence records as a coda::mcmc.list object
setMethod(
          'conv.rec',
          signature=signature(object='pmcmcList'),
          definition=function (object, ...) {
            f <- selectMethod("conv.rec","pmcmc")
            coda::mcmc.list(lapply(object,f,...))
          }
          )

## extract the filtered trajectories from a pmcmc
setMethod(
          'filter.traj',
          signature=signature(object='pmcmc'),
          definition=function (object, ...) {
            f <- selectMethod("filter.traj","pfilterd.pomp")
            f(object,...)
          }
          )

## extract the filtered trajectories from a pmcmcList
setMethod(
          'filter.traj',
          signature=signature(object='pmcmcList'),
          definition=function (object, ...) {
            f <- selectMethod("filter.traj","pmcmc")
            lapply(object,f,...)
          }
          )

## plot pmcmc object
setMethod(
          "plot",
          signature=signature(x='pmcmc'),
          function (x, y, ...) {
            if (!missing(y)) {
              y <- substitute(y)
              warning(sQuote(y)," is ignored")
            }
            pmcmc.diagnostics(list(x))
          }
          )


setMethod(
          "plot",
          signature=signature(x='pmcmcList'),
          definition=function (x, y, ...) {
            if (!missing(y)) {
              y <- substitute(y)
              warning(sQuote(y)," is ignored")
            }
            pmcmc.diagnostics(x)
          }
          )

pmcmc.diagnostics <- function (z) {
  ## assumes that x is a list of pmcmcs with identical structure
  mar.multi <- c(0,5.1,0,2.1)
  oma.multi <- c(6,0,5,0)
  xx <- z[[1]]
  estnames <- xx@pars

  ## plot pmcmc convergence diagnostics
  other.diagnostics <- c("loglik", "log.prior","nfail")
  plotnames <- c(other.diagnostics,estnames)
  nplots <- length(plotnames)
  n.per.page <- min(nplots,10)
  nc <- if (n.per.page<=4) 1 else 2
  nr <- ceiling(n.per.page/nc)
  oldpar <- par(mar=mar.multi,oma=oma.multi,mfcol=c(nr,nc))
  on.exit(par(oldpar)) 
  low <- 1
  hi <- 0
  iteration <- seq(0,xx@Nmcmc)
  while (hi<nplots) {
    hi <- min(low+n.per.page-1,nplots)
    for (i in seq(from=low,to=hi,by=1)) {
      n <- i-low+1
      dat <- sapply(z,conv.rec,pars=plotnames[i])
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
      if (do.xax) mtext("PMCMC iteration",side=1,line=3)
    }  
    low <- hi+1
    mtext("PMCMC convergence diagnostics",3,line=2,outer=TRUE)
  }
  invisible(NULL)
}
