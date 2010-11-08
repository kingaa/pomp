\name{traj.match}
\alias{traj.match}
\alias{traj.match-pomp}
\alias{traj.match,pomp-method}
\alias{traj.match-traj.matched.pomp}
\alias{traj.match,traj.matched.pomp-method}
\alias{logLik,traj.matched.pomp-method}
\alias{logLik-traj.matched.pomp}
\alias{$,traj.matched.pomp-method}
\alias{$-traj.matched.pomp}
\alias{summary,traj.matched.pomp-method}
\alias{summary-traj.matched.pomp}
\alias{traj.matched.pomp-class}
\title{Trajectory matching}
\description{
  Match trajectories to data.
}
\usage{
  \S4method{traj.match}{pomp}(object, start, est,
                 method = c("Nelder-Mead", "SANN", "subplex"),
                 gr = NULL, eval.only = FALSE, \dots)
  \S4method{traj.match}{traj.matched.pomp}(object, start, est,
                 method = c("Nelder-Mead", "SANN", "subplex"),
                 gr = NULL, eval.only = FALSE, \dots)
}
\arguments{
  \item{object}{A \code{pomp} object.}
  \item{start}{initial guess for parameters.}
  \item{est}{character vector containing the names of parameters to be estimated.}
  \item{method}{
    Optimization method.
    Choices are \code{\link[subplex]{subplex}} and any of the methods used by \code{\link{optim}}.
  }
  \item{gr}{
    Passed to \code{\link{optim}}.
  }
  \item{eval.only}{
    logical;
    if \code{TRUE}, no optimization is attempted and the log likelihood value is evaluated at the \code{start} parameters.
  }
  \item{\dots}{
    Arguments that will be passed to \code{\link{optim}} in its \code{control} list.
  }
}
\details{
  Trajectory matching is accomplished using \code{\link{optim}}.
  The \code{\link{trajectory}} method is used for this, which in turn uses the \code{skeleton} slot of the \code{pomp} object.
  The quantity maximized is the likelihood of the data given the trajectory, as returned by \code{\link{dmeasure}}.
}
\value{
  An object of class \code{traj.matched.pomp}.
  This class inherits from class \code{\link[=pomp-class]{pomp}} and contains the following additional slots:
  \describe{
    \item{evals}{
      number of function and gradient evaluations by the optimizer.
      See \code{\link{optim}}.
    }
    \item{value}{
      value of the objective function.
      Larger values indicate better fit (i.e., \code{traj.match} attempts to maximize this quantity.
    }
    \item{convergence, msg}{
      convergence code and message from the optimizer.
      See \code{\link{optim}}.
    }
  }

  Available methods for objects of this type include \code{summary} and \code{logLik}.
  The other slots of this object can be accessed via the \code{$} operator.
}
\examples{
  data(ou2)
  true.p <- c(
	      alpha.1=0.9,alpha.2=0,alpha.3=-0.4,alpha.4=0.99,
	      sigma.1=2,sigma.2=0.1,sigma.3=2,
	      tau=1,
              x1.0=50,x2.0=-50
	      )
  simdata <- simulate(ou2,nsim=1,params=true.p,seed=43553)
  guess.p <- true.p
  res <- traj.match(
		    simdata,
		    start=guess.p,
		    est=c('alpha.1','alpha.3','alpha.4','x1.0','x2.0','tau'),
		    maxit=2000,
		    method="Nelder-Mead",
		    reltol=1e-8
		    )

  summary(res)

  plot(range(time(res)),range(c(obs(res),states(res))),type='n',xlab="time",ylab="x,y")
  points(y1~time,data=as(res,"data.frame"),col='blue')
  points(y2~time,data=as(res,"data.frame"),col='red')
  lines(x1~time,data=as(res,"data.frame"),col='blue')
  lines(x2~time,data=as(res,"data.frame"),col='red')
}
\seealso{
  \code{\link{trajectory}},
  \code{\link{pomp}},
  \code{\link{optim}},
  \code{\link[subplex]{subplex}}
}
\keyword{models}
\keyword{ts}
