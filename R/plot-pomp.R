## plot pomp object

setMethod(
          "plot",
          "pomp",
          function (x, y = NULL, 
                    variables,
                    panel = lines,
                    nc = NULL,
                    yax.flip = FALSE,
                    mar = c(0, 5.1, 0, if(yax.flip) 5.1 else 2.1),
                    oma = c(6, 0, 5, 0),
                    axes = TRUE,
                    ...) {
            X <- as(x,'data.frame')
            vars <- names(X)
            tpos <- match("time",vars)
            if (is.na(tpos))
              stop("'pomp' plot error: no data variable labeled 'time'",call.=FALSE)
            if (missing(variables))
              vars <- vars[-tpos]
            else
              vars <- variables
            plotpomp <- function (x, time, 
                                  xy.labels, xy.lines, panel = lines, nc, xlabel,
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
              if (nser > 10)
                stop("'pomp' plot error: cannot plot more than 10 series as \"multiple\"",call.=FALSE)
              if (is.null(main))
                main <- xlabel
              nm <- colnames(x)
              if (is.null(nm))
                nm <- paste("Series", 1:nser)
              if (is.null(nc))
                nc <- if(nser>4){2}else{1}
              nr <- ceiling(nser/nc)
              oldpar <- par(mar=mar,oma=oma,mfcol=c(nr,nc))
              on.exit(par(oldpar))
              for (i in 1:nser) {
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
            v2 <- min(v1+plots.per.page,length(vars))
            for (page in 1:n.page) {
              vv <- vars[seq(from=v1,to=v2)]
              plotpomp(
                       x=X[vv],
                       time=X[[tpos]],
                       xy.labels=FALSE,
                       xlabel=deparse(substitute(x,env=parent.frame(1))),
                       panel=panel,
                       nc=nc,
                       axes=axes,
                       ...
                       )
              v1 <- v2+1
              v2 <- min(v2+plots.per.page,length(vars))    
            }
            invisible(NULL)
          }
          )
