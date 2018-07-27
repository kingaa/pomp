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
  "values",
  signature=signature(object="probed.pomp"),
  definition=function (object, ...) {
    x <- as.data.frame(rbind(object@datvals,object@simvals))
    row.names(x) <- seq.int(from=0,to=nrow(x)-1)
    x$.id <- factor(c("data",rep("sim",nrow(x)-1)))
    x
  }
)

setMethod("plot",
  signature=signature(x="probed.pomp"),
  definition=function (x, y, ...) {
    if (!missing(y))
      warning("in ",sQuote("plot-probed.pomp"),": ",
        sQuote("y")," is ignored",call.=FALSE)
    probeplot.internal(x=x,...)
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
