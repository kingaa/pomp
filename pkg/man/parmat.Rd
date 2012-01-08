\name{parmat}
\alias{parmat}
\title{Create a matrix of parameters}
\description{
  \code{parmat} is a utility that makes a vector of parameters suitable for use in \pkg{pomp} functions.
}
\usage{
parmat(params, nrep = 1)
}
\arguments{
  \item{params}{
    named numeric vector of parameters.
  }
  \item{nrep}{
    number of replicates (columns) desired.
  }
}
\value{
  \code{parmat} returns a matrix consisting of \code{nrep} copies of \code{params} and rownames identical to the names of \code{params}.
}
\author{Aaron A. King \email{kingaa at umich dot edu}}
\examples{
  ## generate a bifurcation diagram for the Ricker map
  data(ricker)
  p <- parmat(coef(ricker),nrep=500)
  p["log.r",] <- seq(from=1.5,to=4,length=500)
  x <- trajectory(ricker,times=seq(from=1000,to=2000,by=1),params=p)
  matplot(p["log.r",],x["N",,],pch='.',col='black',xlab="log(r)",ylab="N")
}
