compare.mif <- function (z) {
  ## assumes that x is a list of mifs with identical structure
  if (!is.list(z)) z <- list(z)
  if (!all(sapply(z,function(x)is(x,'mif'))))
    stop("compare.mif error: 'z' must be a mif object or a list of mif objects",call.=FALSE)
  mar.multi <- c(0,5.1,0,2.1)
  oma.multi <- c(6,0,5,0)
  xx <- z[[1]]
  ivpnames <- xx@ivps
  estnames <- c(xx@pars,ivpnames)
  parnames <- names(coef(xx))
  unestnames <- parnames[-match(estnames,parnames)]

  ## plot filter means
  filt.diag <- rbind("eff. sample size"=xx@eff.sample.size,filter.mean(xx))
  filtnames <- rownames(filt.diag)
  plotnames <- filtnames[-match(unestnames,filtnames)]
  lognames <- filtnames[1] # eff. sample size
  nplots <- length(plotnames)
  n.per.page <- min(nplots,10)
  if(n.per.page<=4) nc <- 1 else nc <- 2
  nr <- ceiling(n.per.page/nc)
  oldpar <- par(mar=mar.multi,oma=oma.multi,mfcol=c(nr,nc),ask=dev.interactive(orNone=TRUE))
  on.exit(par(oldpar)) 
  low <- 1
  hi <- 0
  time <- time(xx)
  while (hi<nplots) {
    hi <- min(low+n.per.page-1,nplots)
    for (i in low:hi) {
      n <- i-low+1
      logplot <- if (plotnames[i]%in%lognames) "y" else ""
      dat <- sapply(
                    z,
                    function(po, label) {
                      if (label=="eff. sample size") 
                        po@eff.sample.size
                      else
                        filter.mean(po,label)
                    },
                    label=plotnames[i]
                    )
      matplot(
              y=dat, 
              x=time,
              axes = FALSE,
              xlab = "",
              log=logplot,
              ylab = "",
              type = "l"
              )
      box()
      y.side <- 2
      axis(y.side, xpd = NA)
      mtext(plotnames[i], y.side, line = 3)
      do.xax <- (n%%nr==0||n==n.per.page)
      if (do.xax) axis(1,xpd=NA)
      if (do.xax) mtext("time",side=1,line=3)
    }  
    low <- hi+1
    mtext("Filter diagnostics (last iteration)",3,line=2,outer=TRUE)
  } 

  ## plot mif convergence diagnostics
  other.diagnostics <- c("loglik", "nfail")
  plotnames <- c(other.diagnostics,estnames,ivpnames)
  nplots <- length(plotnames)
  n.per.page <- min(nplots,10)
  nc <- if (n.per.page<=4) 1 else 2
  nr <- ceiling(n.per.page/nc)
  par(mar=mar.multi,oma=oma.multi,mfcol=c(nr,nc))
  ## on.exit(par(oldpar)) 
  low <- 1
  hi <- 0
  iteration <- seq(0,xx@Nmif)
  while (hi<nplots) {
    hi <- min(low+n.per.page-1,nplots)
    for (i in low:hi) {
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
      if (do.xax) mtext("MIF iteration",side=1,line=3)
    }  
    low <- hi+1
    mtext("MIF convergence diagnostics",3,line=2,outer=TRUE)
  }
  invisible(NULL)
}


