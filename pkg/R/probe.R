setClass(
         "probed.pomp",
         contains="pomp",
         representation(
                        probes="list",
                        seed="integer",
                        datvals="numeric",
                        simvals="array",
                        quantiles="numeric",
                        pvals="numeric"
                        )
         )

setMethod(
          "summary",
          "probed.pomp",
          function (object, ...) {
            list(
                 coef=coef(object),
                 nsim=nrow(object@simvals),
                 quantiles=object@quantiles,
                 pvals=object@pvals
                 )
          }
          )

setGeneric("probe",function(object,probes,...)standardGeneric("probe"))

setMethod(
          "probe",
          signature(object="pomp"),
          function (object, probes, nsim = 1, seed = NULL, ...) {
            if (!is.list(probes)) probes <- list(probes)
            if (!all(sapply(probes,is.function)))
              stop(sQuote("probes")," must be a function or a list of functions")
            if (!all(sapply(probes,function(f)length(formals(f))==1)))
              stop("each probe must be a function of a single argument")
            if (is.null(seed)) {
              if (exists('.Random.seed',where=.GlobalEnv)) {
                seed <- get(".Random.seed",pos=.GlobalEnv)
              }
            }
            
            ## apply probes to data
            datval <- .Call(apply_probe_data,object,probes)
            ## apply probes to model simulations
            simval <- .Call(
                            apply_probe_sim,
                            object,
                            nsim,
                            coef(object),
                            seed,
                            probes,
                            datval
                            )
                            
            nprobes <- length(datval)
            pvals <- numeric(nprobes)
            names(pvals) <- names(datval)
            quants <- numeric(nprobes)
            names(quants) <- names(datval)
            for (k in seq_len(nprobes)) {
              r <- min(sum(simval[,k]>datval[k]),sum(simval[,k]<datval[k]))
              tails <- (r+1)/(nsim+1)
              pvals[k] <- min(2*tails,1)
              quants[k] <- sum(simval[,k]<datval[k])/nsim
            }

            new(
                "probed.pomp",
                object,
                probes=probes,
                seed=as.integer(seed),
                datvals=datval,
                simvals=simval,
                quantiles=quants,
                pvals=pvals
                )
          }
          )

setMethod(
          "probe",
          signature(object="probed.pomp"),
          function (object, probes, nsim, seed = NULL, ...) {
            if (missing(probes)) probes <- object@probes
            if (!is.list(probes)) probes <- list(probes)
            if (!all(sapply(probes,is.function)))
              stop(sQuote("probes")," must be a function or a list of functions")
            if (is.null(seed)) {
              if (exists('.Random.seed',where=.GlobalEnv)) {
                seed <- get(".Random.seed",pos=.GlobalEnv)
              }
            }
            if (missing(nsim)) nsim <- nrow(object@simvals)
            probe(
                  as(object,"pomp"),
                  probes=probes,
                  nsim=nsim,
                  seed=seed,
                  ...
                  )
          }
          )

setMethod(
          "plot",
          "probed.pomp", 
          function (x, y, ...) {

            ##function for plotting diagonal panels
            diag.panel.hist <- function(x, ...) {
              ##plot a histogram for the simulations
              usr <- par("usr")
              on.exit(par(usr))
              par(usr=c(usr[1:2],0,1.5))
              h <- hist(x[-1],plot=FALSE)
              breaks <- h$breaks
              nB <- length(breaks)
              y <- h$counts
              y <- y/max(y)
              rect(breaks[-nB],0,breaks[-1],y,...)
              ##plot the data point
              lines(c(x[1],x[1]),c(0,max(h$counts)),col="red")
            }

            ##function for plotting above-diagonal panels
            above.diag.panel <- function (x, y, ...) {
              ##plot the simulations
              points(x[-1],y[-1],...)
              ##plot the data
              mMx <- c(min(x),max(x))
              mMy <- c(min(y),max(y))
              lines(c(x[1],x[1]),mMy,col="red")
              lines(mMx,c(y[1],y[1]),col="red")
            }
            
            ##function for plotting below-diagonal panels
            below.diag.panel <- function (x, y, ...) {
              mMx <- c(min(x),max(x))
              mMy <- c(min(y),max(y))
              x <- x[-1]
              y <- y[-1]
              correls <- round(cor(x,y),3)
              text(mean(mMx),mean(mMy),correls,cex=1)
            }
            
            ##prepare the arguments for these functions
            nprobes <- length(x@datvals)
            nsim <- nrow(x@simvals)
            datsimvals <- array(dim=c(nsim+1,nprobes))
            datsimvals[1,] <- x@datvals
            datsimvals[-1,] <- x@simvals
            
            labels <- paste("pb",seq_len(nprobes))
            if (!is.null(names(x@datvals)))
              labels <- ifelse(names(x@datvals)=="",labels,names(x@datvals))
            lab.plus <- paste(labels,paste("p=",round(x@pvals,3),sep=""),sep="\n")
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
          )

