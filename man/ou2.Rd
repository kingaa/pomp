\name{ou2}
\alias{ou2}
\docType{data}
\title{Two-dimensional Ornstein-Uhlenbeck process}
\description{
  \code{ou2} is a \code{pomp} object encoding a 2-D Ornstein-Uhlenbeck process.
}
\usage{data(ou2)}
\details{
  The OU process is fully but noisily observed.
  The functions \code{rprocess}, \code{dprocess}, \code{rmeasure}, and \code{dmeasure} are implemented using compiled C code for computational speed:
  see the source code for details.
  The use of this object is demonstrated in the vignette.
}
\examples{
data(ou2)
plot(ou2)
slotNames(ou2)
p <- c(
       alpha.1=0.9,alpha.2=0,alpha.3=0,alpha.4=0.99,
       sigma.1=1,sigma.2=0,sigma.3=2,
       tau=1,
       x1.0=50,x2.0=-50
       )
x <- simulate(ou2,nsim=10,seed=20348585,params=p)
plot(x[[1]])
pfilter(ou2,params=p,Np=1000)$loglik
}
\seealso{\code{\link{pomp-class}} and the vignettes}
\keyword{datasets}
