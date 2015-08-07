## helper function to extract a covariance matrix from MCMC traces

covmat.internal <- function (traces, start, thin, expand = 2.38, ...) {
  dd <- dim(traces)
  nms <- colnames(traces)
  keep <- seq.int(from=as.integer(start),to=dd[1],by=as.integer(thin))
  if (length(keep) < 100)
    warning("only ",length(keep)," points being used to estimate covariance matrix")
  if (length(dd)==2L) {
    traces <- traces[keep,,drop=FALSE]
  } else if (length(dd)==3L) {
    traces <- aperm(traces[keep,,,drop=FALSE],c(1L,3L,2L))
    dd <- dim(traces)
    dim(traces) <- c(dd[1L]*dd[2L],dd[3L])
  }
  v <- var(traces)
  keep <- which(diag(v)>0)
  v <- expand^2*v[keep,keep]/length(keep)
  dimnames(v) <- list(nms[keep],nms[keep])
  v
}

setMethod("covmat",signature=signature(object="pmcmc"),
          definition=function (object, start = 1, thin = 1,
            expand = 2.38, ...) {
            covmat.internal(traces=as.matrix(conv.rec(object,object@pars)),
                            start=start,thin=thin,
                            expand=expand)
          })

setMethod("covmat",signature=signature(object="pmcmcList"),
          definition=function (object, start = 1, thin = 1,
            expand = 2.38, ...) {
            pars <- do.call(unique,lapply(object,slot,"pars"))
            covmat.internal(traces=as.array(conv.rec(object,pars)),
                            start=start,thin=thin,
                            expand=expand)
          })

setMethod("covmat",signature=signature(object="abc"),
          definition=function (object, start = 1, thin = 1,
            expand = 2.38, ...) {
            covmat.internal(traces=as.matrix(conv.rec(object,object@pars)),
                            start=start,thin=thin,
                            expand=expand)
          })

setMethod("covmat",signature=signature(object="abcList"),
          definition=function (object, start = 1, thin = 1,
            expand = 2.38, ...) {
            pars <- do.call(unique,lapply(object,slot,"pars"))
            covmat.internal(traces=as.array(conv.rec(object,pars)),
                            start=start,thin=thin,
                            expand=expand)
          })
