slice.design <- function (vars, n) {
  if (!is.list(vars) || is.null(names(vars)))
    stop("slice.design error: ",sQuote("vars")," must be a named list")
  if (!all(sapply(vars,function(x)is.numeric(x)&&(length(x)%in%c(1,3)))))
    stop("slice.design error: each entry in ",sQuote("vars")," must specify a center or a center and range")
  nvars <- length(vars)
  varying <- which(sapply(vars,length)>1)
  nranges <- length(varying)
  center.point <- sapply(vars,function(x)if(length(x)>1){x[2]}else{x[1]})
  ranges <- lapply(vars[varying],function(x)x[c(1,3)])
  x <- as.data.frame(matrix(center.point,byrow=TRUE,ncol=nvars,nrow=n,dimnames=list(NULL,names(vars))))
  y <- vector(mode="list",length=nranges)
  for (v in seq(length=nranges)) {
    y[[v]] <- x
    w <- varying[v]
    y[[v]][[w]] <- seq(from=ranges[[v]][1],to=ranges[[v]][2],length=n)
  }
  y <- do.call(rbind,y)
  y$slice <- as.factor(rep(names(vars)[varying],each=n))
  y
}

profile.design <- function (..., vars, n) {
  prof <- list(...)
  y <- as.matrix(sobol(vars,n))
  x <- as.matrix(do.call(expand.grid,prof))
  z <- array(dim=c(nrow(x),nrow(y),ncol(x)+ncol(y)))
  for (j in 1:nrow(x)) {
    for (k in 1:nrow(y)) {
      z[j,k,] <- c(x[j,],y[k,])
    }
  }
  dim(z) <- c(nrow(x)*nrow(y),ncol(x)+ncol(y))
  colnames(z) <- c(colnames(x),colnames(y))
  as.data.frame(z)
}
