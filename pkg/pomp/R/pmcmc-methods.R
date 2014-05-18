## this file contains short definitions of methods for the 'pmcmc' class

## extract the estimated log likelihood
setMethod('logLik','pmcmc',function(object,...)object@loglik)

## extract the filtering means
setMethod(
          'filter.mean',
          'pmcmc',
          function (object, pars, ...) {
            if (missing(pars)) pars <- rownames(object@filter.mean)
            object@filter.mean[pars,]
          }
          )

## extract the convergence record
setMethod(
          'conv.rec',
          'pmcmc',
          function (object, pars, ...) {
            if (missing(pars)) pars <- colnames(object@conv.rec)
            object@conv.rec[,pars]
          }
          )

## plot pmcmc object
setMethod(
          "plot",
          "pmcmc",
          function (x, y = NULL, ...) {
            compare.pmcmc(x)
          }
          )

compare.pmcmc <- function (z) {
  ## assumes that x is a list of pmcmcs with identical structure
  if (!is.list(z)) z <- list(z)
  if (!all(sapply(z,function(x)is(x,'pmcmc'))))
    stop("compare.pmcmc error: ",sQuote("z"),
         " must be a pmcmc object or a list of pmcmc objects",call.=FALSE)
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
      if (do.xax) mtext("PMCMC iteration",side=1,line=3)
    }  
    low <- hi+1
    mtext("PMCMC convergence diagnostics",3,line=2,outer=TRUE)
  }
  invisible(NULL)
}
