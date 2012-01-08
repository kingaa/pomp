parmat <- function (params, nrep = 1) {
  p <- matrix(data=params,nrow=length(params),ncol=nrep)
  rownames(p) <- names(params)
  p
}
