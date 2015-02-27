\name{verhulst}
\alias{verhulst}
\docType{data}
\title{Simple Verhulst-Pearl (logistic) model.}
\description{
  \code{verhulst} is a \code{pomp} object encoding a univariate stochastic logistic model with measurement error.
}
\details{
  The model is written as an Ito diffusion, \eqn{dn = r n (1-n/K) dt + \sigma n dW}, where \eqn{W} is a Wiener process.
  It is implemented using the \code{\link{euler.sim}} plug-in.
}
\examples{
pompExample(verhulst)
plot(verhulst)
coef(verhulst)
params <- cbind(
                c(n.0=100,K=10000,r=0.2,sigma=0.4,tau=0.1),
                c(n.0=1000,K=11000,r=0.1,sigma=0.4,tau=0.1)
                )
x <- simulate(verhulst,params=params,states=TRUE)
matplot(time(verhulst),t(x['n',,]),type='l')
y <- trajectory(verhulst,params=params)
matlines(time(verhulst),t(y['n',,]),type='l',lwd=2)
}
\seealso{
  \code{\link{pomp}} and the tutorial documents on the \href{http://pomp.r-forge.r-project.org}{package website}.
}
\keyword{models}
\keyword{datasets}
