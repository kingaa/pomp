##' pomp plotting facilities
##'
##' Diagnostic plots.
##'
##' @name plot
##' @rdname plot
##' @include pomp_class.R
##' @include abc.R mif2.R pmcmc.R pfilter.R spect.R probe.R wpfilter.R
##' @include listie.R
##' @aliases plot,missing-method
##' @importFrom graphics par abline pairs matplot box axis mtext points polygon lines plot.default legend hist rect text title
##' @importFrom grDevices rgb dev.interactive
##' @importFrom stats quantile cor density
##'
NULL

setGeneric("plot")

setClassUnion("pomp_plottable",c("pomp","pfilterd_pomp","wpfilterd_pomp"))

##' @rdname plot
##' @param x the object to plot
##' @param variables optional character; names of variables to be displayed
##' @param panel function of prototype \code{panel(x, col, bg, pch, type, ...)} which gives the action to be carried out in each panel of the display.
##' @param nc the number of columns to use.
##' Defaults to 1 for up to 4 series, otherwise to 2.
##' @param yax.flip logical;
##' if TRUE, the y-axis (ticks and numbering) should flip from side 2 (left) to 4 (right) from series to series.
##' @param mar,oma the \code{\link{par}} \code{mar} and \code{oma} settings.
##' Modify with care!
##' @param axes logical; indicates if x- and y- axes should be drawn
##' @param \dots ignored or passed to low-level plotting functions
##' @export
setMethod(
  "plot",
  signature=signature(x="pomp_plottable"),
  definition=function (x, variables,
    panel = lines, nc = NULL, yax.flip = FALSE,
    mar = c(0, 5.1, 0, if (yax.flip) 5.1 else 2.1),
    oma = c(6, 0, 5, 0), axes = TRUE, ...) {

    plotpomp.internal(x=x,variables=variables,
      panel=panel,nc=nc,yax.flip=yax.flip,
      mar=mar,oma=oma,axes=axes,...)

  }
)

##' @rdname plot
##' @param pars names of parameters.
##' @export
setMethod(
  "plot",
  signature=signature(x="Pmcmc"),
  definition=function (x, ..., pars) {
    plot(traces(x,pars),...)
  }
)

##' @rdname plot
##' @param scatter logical; if \code{FALSE}, traces of the parameters named in \code{pars} will be plotted against ABC iteration number.
##' If \code{TRUE}, the traces will be displayed or as a scatterplot.
##' @export
setMethod(
  "plot",
  signature=signature(x="Abc"),
  definition=function (x, ..., pars, scatter = FALSE) {
    abc.diagnostics(x,pars=pars,scatter=scatter,...)
  }
)

##' @rdname plot
##' @param transform logical; should the parameter be transformed onto the estimation scale?
##' @export
setMethod(
  "plot",
  signature=signature(x="Mif2"),
  definition=function (x, ..., pars, transform = FALSE) {
    mif2.diagnostics(x,pars,transform=as.logical(transform))
  }
)

##' @rdname plot
##' @param y ignored
##' @export
setMethod(
  "plot",
  signature=signature(x="probed_pomp"),
  definition=function (x, ...) {
    probeplot.internal(x,...)
  }
)

##' @rdname plot
##' @param max.plots.per.page positive integer; maximum number of plots on a page
##' @param plot.data logical; should the data spectrum be included?
##' @param quantiles numeric; quantiles to display
##' @param quantile.styles list; plot styles to use for quantiles
##' @param data.styles list; plot styles to use for data
##' @export
setMethod(
  "plot",
  signature=signature(x="spectd_pomp"),
  definition=function (x, ..., max.plots.per.page = 4, plot.data = TRUE,
    quantiles = c(.025, .25, .5, .75, .975),
    quantile.styles = list(lwd=1, lty=1, col="gray70"),
    data.styles = list(lwd=2, lty=2, col="black")) {

    tryCatch(
      plot_spect.internal(x,max.plots.per.page=max.plots.per.page,
        plot.data=plot.data,quantiles=quantiles,quantile.styles=quantile.styles,
        data.styles=data.styles),
      error = function (e) pStop("plot",conditionMessage(e))
    )

  }
)

plotpomp.internal <- function (x, variables,
  panel = lines, nc = NULL, yax.flip = FALSE,
  mar = c(0, 5.1, 0, if (yax.flip) 5.1 else 2.1),
  oma = c(6, 0, 5, 0), axes = TRUE, ...) {

  X <- as(x,"data.frame")
  vars <- names(X)
  tpos <- match(x@timename,vars)
  if (missing(variables)) {
    vars <- vars[-tpos]
    vars <- setdiff(vars,get_covariate_names(x@covar))
    ylabels <- NULL
  } else {
    vars <- variables
    ylabels <- names(variables)
  }

  plotpomp <- function (x, time,
    xy.labels, xy.lines, panel = lines, nc,
    type = "l", xlim = NULL, ylim = NULL, xlab = "time",
    ylab, log = "", col = par("col"), bg = NA,
    pch = par("pch"),
    cex = par("cex"), lty = par("lty"), lwd = par("lwd"),
    axes = TRUE, frame.plot = axes, ann = par("ann"),
    main = NULL,
    ...) {
    panel <- match.fun(panel)
    addmain <- function (main,
      cex.main = par("cex.main"),
      font.main = par("font.main"),
      col.main = par("col.main"),
      ...) {
      mtext(main,side=3,line=3,cex=cex.main,font=font.main,col=col.main,...)
    }
    nser <- NCOL(x)
    nm <- ylab
    if (is.null(nm)) nm <- colnames(x)
    if (is.null(nc)) nc <- if(nser>4){2}else{1}
    nr <- ceiling(nser/nc)
    oldpar <- par(mar=mar,oma=oma,mfcol=c(nr,nc))
    on.exit(par(oldpar))
    for (i in seq_len(nser)) {
      plot.default(y=x[[i]],x=time,axes=FALSE,xlab="",ylab="",log=log,
        col=col,bg=bg,pch=pch,ann=ann,type="n",...)
      panel(y=x[[i]],x=time,col=col,bg=bg,pch=pch,type=type,...)
      if (frame.plot)
        box(...)
      y.side <- if(i%%2||!yax.flip){2}else{4}
      do.xax <- (i%%nr==0)||(i==nser)
      if (axes) {
        axis(y.side, xpd = NA)
        if (do.xax)
          axis(1, xpd = NA)
      }
      if (ann) {
        mtext(nm[i], y.side, line = 3, ...)
        if (do.xax)
          mtext(xlab, side = 1, line = 3, ...)
      }
    }
    if (ann && !is.null(main)) {
      par(mfcol=c(1,1))
      addmain(main,...)
    }
    invisible(NULL)
  }
  n.page <- ceiling(length(vars)/10)
  plots.per.page <- ceiling(length(vars)/n.page)
  if (n.page > 1) {
    op <- par(ask=dev.interactive(orNone=TRUE))
    on.exit(par(op))
  }
  v1 <- 1
  v2 <- min(v1+plots.per.page-1,length(vars))
  for (page in seq_len(n.page)) {
    vv <- vars[seq(from=v1,to=v2)]
    plotpomp(
      x=X[vv],
      time=X[[tpos]],
      xlab=x@timename,
      ylab=ylabels,
      xy.labels=FALSE,
      panel=panel,
      nc=nc,
      axes=axes,
      ...
    )
    v1 <- v2+1
    v2 <- min(v1+plots.per.page-1,length(vars))
  }
  invisible(NULL)
}

abc.diagnostics <- function (z, pars, scatter = FALSE, ...) {
  if (!is.list(z)) z <- list(z)
  if (missing(pars)) {
    pars <- unique(do.call(c,lapply(z,slot,"pars")))
    if (length(pars)<1)
      pars <- unique(do.call(c,lapply(z,function(x)names(x@params))))
  }
  if (scatter) {
    x <- lapply(z,function(x)as.matrix(traces(x,pars)))
    x <- lapply(seq_along(x),function(n)cbind(x[[n]],.num=n))
    x <- do.call(rbind,x)
    if (ncol(x)<3) {
      pStop("plot","can't make a scatterplot with only one variable.")
    } else {
      pairs(x[,pars],col=x[,".num"],...)
    }
  } else {
    mar.multi <- c(0,5.1,0,2.1)
    oma.multi <- c(6,0,5,0)
    xx <- z[[1]]
    estnames <- pars
    ## plot abc convergence diagnostics
    other.diagnostics <- c()
    plotnames <- c(other.diagnostics,estnames)
    nplots <- length(plotnames)
    n.per.page <- min(nplots,10)
    nc <- if (n.per.page<=4) 1 else 2
    nr <- ceiling(n.per.page/nc)
    oldpar <- par(mar=mar.multi,oma=oma.multi,mfcol=c(nr,nc))
    on.exit(par(oldpar))
    low <- 1
    hi <- 0
    iteration <- seq(0,xx@Nabc)
    while (hi<nplots) {
      hi <- min(low+n.per.page-1,nplots)
      for (i in seq(from=low,to=hi,by=1)) {
        n <- i-low+1
        dat <- sapply(z,traces,pars=plotnames[i])
        matplot(y=dat,x=iteration,axes = FALSE,xlab = "",ylab = "",type = "l")
        box()
        y.side <- 2
        axis(y.side,xpd=NA)
        mtext(plotnames[i],y.side,line=3)
        do.xax <- (n%%nr==0||n==n.per.page)
        if (do.xax) axis(1,xpd=NA)
        if (do.xax) mtext("ABC iteration",side=1,line=3)
      }
      low <- hi+1
      mtext("ABC convergence diagnostics",3,line=2,outer=TRUE)
    }
  }
  invisible(NULL)
}

mif2.diagnostics <- function (z, pars, transform) {
  if (!is.list(z)) z <- list(z)
  ## assumes that z is a list of mif2d_pomps with identical structure
  mar.multi <- c(0,5.1,0,2.1)
  oma.multi <- c(6,0,5,0)
  xx <- z[[1]]
  estnames <- names(coef(xx,pars=pars,transform=transform))

  ## plot ESS and cond.logLik
  nplots <- 2
  n.per.page <- 2
  nc <- 1
  nr <- 2
  oldpar <- par(mar=mar.multi,oma=oma.multi,mfcol=c(2,1),
    ask=dev.interactive(orNone=TRUE))
  on.exit(par(oldpar))

  time <- time(xx)
  dat <- sapply(z,eff.sample.size)
  matplot(y=dat,x=time,axes = FALSE,xlab = "",log="y",
    ylab = "eff. sample size",type = "l")
  box()
  axis(2, xpd = NA)
  mtext("Filter diagnostics (last iteration)",side=3,line=2,outer=TRUE)

  dat <- sapply(z,cond.logLik)
  matplot(
    y=dat,
    x=time,
    axes = FALSE,
    xlab = "",
    log='',
    ylab = "cond. logLik",
    type = "l"
  )
  box()
  axis(2, xpd = NA)

  axis(1,xpd=NA)
  mtext("time",side=1,line=3)

  ## plot mif convergence diagnostics
  other.diagnostics <- c("loglik")
  plotnames <- c(other.diagnostics,estnames)
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
    for (i in seq(from=low,to=hi,by=1)) {
      n <- i-low+1
      matplot(
        y=sapply(
          z,
          function (po, label) {
            traces(
              po,label,
              transform=(transform && label %in% estnames)
            )
          },
          label=plotnames[i]
        ),
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
    mtext("MIF2 convergence diagnostics",3,line=2,outer=TRUE)
  }
  invisible(NULL)
}

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

  labels <- paste("probe",seq_len(nprobes))
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

plot_spect.internal <- function (x, max.plots.per.page, plot.data,
  quantiles, quantile.styles, data.styles) {

  spomp <- x
  nquants <- length(quantiles)

  if (!is.list(quantile.styles))
    pStop_(sQuote("quantile.styles")," must be a list.")

  for (i in c("lwd", "lty", "col")) {
    if (is.null(quantile.styles[[i]]))
      quantile.styles[[i]] <- rep(1,nquants)
    if (length(quantile.styles[[i]])==1)
      quantile.styles[[i]] <- rep(quantile.styles[[i]],nquants)
    if (length(quantile.styles[[i]])<nquants) {
      pWarn("plot",sQuote("quantile.styles"),
        " contains an element with more than 1 entry but ",
        "fewer entries than quantiles.")
      quantile.styles[[i]]<-rep(quantile.styles[[i]],nquants)
    }
  }

  if (plot.data) {
    nreps <- ncol(spomp@datspec)

    if (!is.list(data.styles))
      pStop_(sQuote("data.styles")," must be a list")

    for (i in c("lwd", "lty", "col")) {
      if(is.null(data.styles[[i]]))
        data.styles[[i]] <- rep(2,nreps)
      if(length(data.styles[[i]])==1)
        data.styles[[i]] <- rep(data.styles[[i]],nreps)
      if(length(data.styles[[i]]) < nreps) {
        pWarn("plot",sQuote("data.styles"),
          " contains an element with more than 1 entry but ",
          "fewer entries than observed variables.")
        data.styles[[i]] <- rep(data.styles[[i]],nreps)
      }
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
