parmat <- function (params, nrep = 1) {
  params <- as.matrix(params)
  p <- matrix(data=params,nrow=nrow(params),ncol=ncol(params)*nrep)
  rownames(p) <- rownames(params)
  p
}
