\name{traj.match}
\alias{traj.match}
\title{Trajectory matching}
\description{
  Match trajectories.
}
\usage{
traj.match(object, start, est, method = "Nelder-Mead", gr = NULL, \dots)
}
\arguments{
  \item{object}{A \code{pomp} object.}
  \item{start}{Guessed parameters.}
  \item{est}{Character vector containing the names of parameters to be estimated.}
  \item{method}{}
  \item{gr}{}
  \item{\dots}{Arguments that will be passed to \code{optim} in the
    \code{control} list.}
}
\details{
  Trajectory matching is accomplished using \code{\link{optim}}.  It is
  assumed that the process model is deterministic.
}
\examples{
  data(ou2)
  true.p <- c(
	      alpha.1=0.9,alpha.2=0,alpha.3=0,alpha.4=0.99,
	      sigma.1=1,sigma.2=0,sigma.3=2,
	      tau=1,
              x1.0=50,x2.0=-50
	      )
  simdata <- simulate(ou2,nsim=1,params=true.p,seed=43553)
  guess.p <- true.p
  guess.p[grep('sigma',names(guess.p))] <- 0
  res <- traj.match(
		    simdata,
		    start=guess.p,
		    est=c('alpha.1','alpha.4','x1.0','x2.0','tau'),
		    maxit=2000,
		    reltol=1e-8
		    )
}
\seealso{
  \code{\link{optim}},
  \code{\link{pomp}},
  \code{\link{pomp-class}}
}
\keyword{models}
\keyword{ts}
