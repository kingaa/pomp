## this file contains short definitions of methods for the 'abc' class

## abcList class
setClass(
  "abcList",
  contains="list",
  validity=function (object) {
    if (length(object) > 0) {
      if (!all(vapply(object,is,logical(1),"abcd_pomp"))) {
        retval <- paste0(
          "error in ",sQuote("c"),
          ": dissimilar objects cannot be combined"
        )
        return(retval)
      }
      d <- sapply(object,function(x)dim(x@traces))
      if (!all(apply(d,1,diff)==0)) {
        retval <- paste0(
          "error in ",sQuote("c"),
          ": to be combined, ",sQuote("abcd_pomp"),
          " objects must have chains of equal length"
        )
        return(retval)
      }
    }
    TRUE
  }
)

setClassUnion("Abc",c("abcd_pomp","abcList"))

setMethod(
  "concat",
  signature=signature(...="Abc"),
  definition=function (...) {
    y <- lapply(
      list(...),
      function (z) {
        if (is(z,"list"))
          setNames(as(z,"list"),names(z))
        else
          z
      }
    )
    new("abcList",unlist(y))
  }
)

c.Abc <- concat

setMethod(
  "[",
  signature=signature(x="abcList"),
  definition=function(x, i, ...) {
    y <- as(x,"list")
    names(y) <- names(x)
    y <- unlist(y[i])
    if (is.null(y)) {
      list(NULL)
    } else {
      new("abcList",y)
    }
  }
)

## extract the convergence record as an 'mcmc' object
setMethod(
  "traces",
  signature=signature(object="abcd_pomp"),
  definition=function (object, pars, ...) {
    if (missing(pars)) pars <- colnames(object@traces)
    coda::mcmc(object@traces[,pars,drop=FALSE])
  }
)

## extract the convergence record as an 'mcmc.list' object
setMethod(
  "traces",
  signature=signature(object="abcList"),
  definition=function (object, ...) {
    coda::mcmc.list(lapply(object,traces,...))
  }
)

setMethod(
  "show",
  signature=signature(object="abcd_pomp"),
  definition=function (object) {
    cat("<object of class ",sQuote("abcd_pomp"),">\n",sep="")
    invisible(NULL)
  }
)

setMethod(
  "show",
  signature=signature(object="abcList"),
  definition=function (object) {
    y <- as(object,"list")
    names(y) <- names(object)
    show(y)
  }
)

setMethod(
  "plot",
  signature=signature(x="Abc"),
  definition=function (x, y, pars, scatter = FALSE, ...) {
    if (!missing(y))
      warning("in ",sQuote("plot"),": ",sQuote("y")," is ignored",call.=FALSE)
    abc.diagnostics(x,pars=pars,scatter=scatter,...)
  }
)

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
      stop("in ",sQuote("plot-abc"),
        ": can't make a scatterplot with only one variable",call.=FALSE)
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
        if (do.xax) mtext("ABC iteration",side=1,line=3)
      }
      low <- hi+1
      mtext("ABC convergence diagnostics",3,line=2,outer=TRUE)
    }
  }
  invisible(NULL)
}
