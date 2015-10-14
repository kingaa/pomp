## this file contains short definitions of methods for the 'abc' class

## abcList class
setClass(
    'abcList',
    contains='list',
    validity=function (object) {
        if (!all(sapply(object,is,'abc'))) {
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
                ": to be combined, ",sQuote("abc"),
                " objects must have chains of equal length"
            )
            return(retval)
        }
        TRUE
    }
)

setMethod(
    'c',
    signature=signature(x='abc'),
    definition=function (x, ...) {
        y <- list(...)
        if (length(y)==0) {
            new("abcList",list(x))
        } else {
            p <- sapply(y,is,'abc')
            pl <- sapply(y,is,'abcList')
            if (!all(p||pl))
                stop("cannot mix ",sQuote("abc"),
                     " and non-",sQuote("abc")," objects")
            y[p] <- lapply(y[p],list)
            y[pl] <- lapply(y[pl],as,"list")
            new("abcList",c(list(x),y,recursive=TRUE))
        }
    }
)

setMethod(
    'c',
    signature=signature(x='abcList'),
    definition=function (x, ...) {
        y <- list(...)
        if (length(y)==0) {
            x
        } else {
            p <- sapply(y,is,'abc')
            pl <- sapply(y,is,'abcList')
            if (!all(p||pl))
                stop("cannot mix ",sQuote("abc"),
                     " and non-",sQuote("abc")," objects")
            y[p] <- lapply(y[p],list)
            y[pl] <- lapply(y[pl],as,"list")
            new("abcList",c(as(x,"list"),y,recursive=TRUE))
        }
    }
)

setMethod(
    "[",
    signature=signature(x="abcList"),
    definition=function(x, i, ...) {
        new('abcList',as(x,"list")[i])
    }
)

## extract the convergence record as an 'mcmc' object
setMethod(
    'conv.rec',
    signature=signature(object='abc'),
    definition=function (object, pars, ...) {
        if (missing(pars)) pars <- colnames(object@conv.rec)
        coda::mcmc(object@conv.rec[,pars,drop=FALSE])
    }
)

## extract the convergence record as an 'mcmc.list' object
setMethod(
    'conv.rec',
    signature=signature(object='abcList'),
    definition=function (object, ...) {
        f <- selectMethod("conv.rec","abc")
        coda::mcmc.list(lapply(object,f,...))
    }
)

## plot abc object
setMethod(
    "plot",
    signature=signature(x="abc"),
    definition=function (x, y, pars, scatter = FALSE, ...) {
        if (!missing(y)) {
            y <- substitute(y)
            warning(sQuote(y)," is ignored")
        }
        abc.diagnostics(c(x),pars=pars,scatter=scatter,...)
    }
)

setMethod(
    "plot",
    signature=signature(x='abcList'),
    definition=function (x, y, ...) {
        if (!missing(y)) {
            y <- substitute(y)
            warning(sQuote(y)," is ignored")
        }
        abc.diagnostics(x,...)
    }
)

abc.diagnostics <- function (z, pars, scatter = FALSE, ...) {
    if (missing(pars))
        pars <- unique(do.call(c,lapply(z,slot,'pars')))

    if (scatter) {

        x <- lapply(z,function(x)as.matrix(conv.rec(x,pars)))
        x <- lapply(seq_along(x),function(n)cbind(x[[n]],.num=n))
        x <- do.call(rbind,x)
        if (ncol(x)<3) {
            stop("can't make a scatterplot with only one variable")
        } else {
            pairs(x[,pars],col=x[,'.num'],...)
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
                if (do.xax) mtext("ABC iteration",side=1,line=3)
            }  
            low <- hi+1
            mtext("ABC convergence diagnostics",3,line=2,outer=TRUE)
        }

    }
    invisible(NULL)
}
