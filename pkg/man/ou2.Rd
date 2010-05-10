\name{ou2}
\alias{ou2}
\docType{data}
\title{Two-dimensional Ornstein-Uhlenbeck process}
\description{
  \code{ou2} is a \code{pomp} object encoding a bivariate discrete-time Ornstein-Uhlenbeck process.
}
\usage{data(ou2)}
\details{
  The process is fully but noisily observed.
  The functions \code{rprocess}, \code{dprocess}, \code{rmeasure}, \code{dmeasure}, and \code{skeleton} are implemented using compiled C code for computational speed:
  see the source code for details.
  The use of this object is demonstrated in the vignette "intro_to_pomp".
}
\examples{
data(ou2)
plot(ou2)
coef(ou2)
theta <- c(
           alpha.1=0.8, alpha.2=-0.5, alpha.3=0.3, alpha.4=0.9,
           sigma.1=3, sigma.2=-0.5, sigma.3=2,
           tau=1, 
           x1.0=-3, x2.0=4
           )
x <- simulate(ou2,nsim=10,seed=20348585,params=theta)
plot(x[[1]])
pfilter(ou2,params=theta,Np=1000)$loglik
}
\seealso{\code{\link{pomp}} and the vignettes}
\keyword{datasets}
