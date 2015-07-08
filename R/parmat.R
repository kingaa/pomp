parmat <- function (params, nrep = 1) {
  params <- as.matrix(params)
  matrix(data=params,nrow=nrow(params),ncol=ncol(params)*nrep,
         dimnames=list(variable=rownames(params),rep=NULL))
}
