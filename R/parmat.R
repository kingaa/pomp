##' Create a matrix of parameters
##'
##' \code{parmat} is a utility that makes a vector of parameters suitable for
##' use in \pkg{pomp} functions.
##'
##'
##' @param params named numeric vector or matrix of parameters.
##' @param nrep number of replicates (columns) desired.
##'
##' @return \code{parmat} returns a matrix consisting of \code{nrep} copies of
##' \code{params}.
##'
##' @author Aaron A. King
##'
##' @examples
##'
##'   ## generate a bifurcation diagram for the Ricker map
##'   p <- parmat(coef(ricker()),nrep=500)
##'   p["r",] <- exp(seq(from=1.5,to=4,length=500))
##'   x <- trajectory(ricker(),times=seq(from=1000,to=2000,by=1),params=p)
##'   matplot(p["r",],x["N",,],pch='.',col='black',xlab="log(r)",ylab="N",log='x')
##'
##' @export

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
