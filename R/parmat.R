parmat <- function (params, nrep = 1) {
  d <- dim(params)
  if (is.null(d) || length(d) == 1) {
    matrix(data=params,nrow=length(params),ncol=nrep,
      dimnames=list(variable=names(params),rep=NULL))
  } else if (length(d) == 2) {
    matrix(data=params,nrow=nrow(params),ncol=ncol(params)*nrep,
      dimnames=list(variable=rownames(params),rep=NULL))
  } else {
    matrix(data=params,nrow=nrow(params),ncol=prod(d[-1])*nrep,
      dimnames=list(variable=rownames(params),rep=NULL))
  }
}
