## this file contains short definitions of methods for the 'abc' class

## extract the convergence record
setMethod(
          'conv.rec',
          'abc',
          function (object, pars, ...) {
            if (missing(pars)) pars <- colnames(object@conv.rec)
            coda::mcmc(object@conv.rec[,pars])
          }
          )

## plot abc object
setMethod(
          "plot",
          "abc",
          function (x, y, pars, scatter = FALSE, ...) {
            if (missing(pars)) pars <- x@pars
            if (scatter) {
              pairs(as.matrix(conv.rec(x, pars)))
            } else {
              plot.ts(conv.rec(x,pars),xlab="iteration",...)
            }
          }
          )

compare.abc <- function (z) {
  ## assumes that z is a list of abcs with identical structure
  if (!is.list(z)) z <- list(z)
  if (!all(sapply(z,function(x)is(x,'abc'))))
    stop("compare.abc error: ",
         sQuote("z"),
         " must be a pmcmc object or a list of pmcmc objects",call.=FALSE)
  mar.multi <- c(0,5.1,0,2.1)
  oma.multi <- c(6,0,5,0)
  xx <- z[[1]]
  estnames <- xx@pars
  parnames <- names(coef(xx))
  unestnames <- parnames[-match(estnames,parnames)]

  ## plot pmcmc convergence diagnostics
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
      if (do.xax) mtext("ABC iteration",side=1,line=3)
    }  
    low <- hi+1
    mtext("ABC convergence diagnostics",3,line=2,outer=TRUE)
  }
  invisible(NULL)
}
