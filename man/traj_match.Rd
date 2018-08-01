\name{Trajectory matching}
\title{Parameter estimation by fitting the trajectory of a model's deterministic skeleton to data}
\alias{Trajectory matching}
\alias{traj.match}
\alias{traj.match-pomp}
\alias{traj.match,ANY-method}
\alias{traj.match,missing-method}
\alias{traj.match,pomp-method}
\alias{traj.match-traj_matched_pomp}
\alias{traj.match,traj_matched_pomp-method}
\alias{logLik,traj_matched_pomp-method}
\alias{logLik-traj_matched_pomp}
\alias{$,traj_matched_pomp-method}
\alias{$-traj_matched_pomp}
\alias{summary,traj_matched_pomp-method}
\alias{summary-traj_matched_pomp}
\alias{traj_matched_pomp-class}
\alias{traj.match.objfun}
\alias{traj.match.objfun-pomp}
\alias{traj.match.objfun,ANY-method}
\alias{traj.match.objfun,missing-method}
\alias{traj.match.objfun,pomp-method}
\description{
  This function attempts to match trajectories of a model's deterministic skeleton to data.
  Trajectory matching is equivalent to maximum likelihood estimation under the assumption that process noise is entirely absent,
  i.e., that all stochasticity is measurement error.
  Accordingly, this method uses only the \code{skeleton} and \code{dmeasure} components of a \acronym{POMP} model.
}
\usage{
  \S4method{traj.match}{pomp}(object, start, est = character(0),
           method = c("Nelder-Mead","subplex","SANN","BFGS",
                      "sannbox","nloptr"),
           transform = FALSE, \dots)
  \S4method{traj.match}{traj_matched_pomp}(object, est, transform, \dots)
  \S4method{traj.match.objfun}{pomp}(object, params, est, transform = FALSE, \dots)
}
\arguments{
  \item{object}{
    A \code{\link{pomp}} object.
    If \code{object} has no \code{skeleton} slot, an error will be generated.
  }
  \item{start}{
    named numeric vector containing an initial guess for parameters.
    By default \code{start=coef(object)} if the latter exists.
  }
  \item{params}{
    optional named numeric vector of parameters.
    This should contain all parameters needed by the \code{skeleton} and \code{dmeasure} slots of \code{object}.
    In particular, any parameters that are to be treated as fixed should be present here.
    Parameter values given in \code{params} for parameters named in \code{est} will be ignored.
    By default, \code{params=coef(object)} if the latter exists.
  }
  \item{est}{
    character vector containing the names of parameters to be estimated.
    In the case of \code{traj.match.objfun}, the objective function that is constructed will assume that its argument contains the parameters in this order.
  }
  \item{method}{
    Optimization method.
    Choices are \code{\link[subplex]{subplex}}, \dQuote{sannbox}, any of the methods used by \code{\link{optim}}, and \code{nloptr}.
    The latter makes available all the optimization algorithms of the \pkg{nloptr} package (\url{https://cran.r-project.org/package=nloptr}).
  }
  \item{transform}{
    logical;
    if \code{TRUE}, optimization is performed on the transformed scale.
  }
  \item{\dots}{
    Extra arguments that will be passed either to the optimizer (\code{\link{optim}}, \code{\link[subplex]{subplex}}, \code{\link[nloptr]{nloptr}}, or \code{\link{sannbox}}, via their \code{control} (\code{optim}, \code{subplex}, \code{sannbox}) or \code{opts} (\code{nloptr}) lists) or to the ODE integrator.
    In \code{traj.match}, extra arguments will be passed to the optimizer.
    In \code{traj.match.objfun}, extra arguments are passed to \code{\link{trajectory}}.
    If extra arguments are needed by both optimizer and \code{\link{trajectory}}, construct an objective function first using \code{traj.match.objfun}, then give this objective function to the optimizer.
  }
}
\details{
  In \pkg{pomp}, trajectory matching is the term used for maximizing the likelihood of the data under the assumption that there is no process noise.
  Specifically, \code{traj.match} calls an optimizer (\code{\link{optim}}, \code{\link[subplex]{subplex}}, and \code{\link{sannbox}} are the currently supported options) to minimize an objective function.
  For any value of the model parameters, this objective function is calculated by
  \enumerate{
    \item computing the deterministic trajectory of the model given the parameters.
    This is the trajectory returned by \code{\link{trajectory}}, which relies on the model's deterministic skeleton as specified in the construction of the \code{pomp} object \code{object}.
    \item evaluating the negative log likelihood of the data under the measurement model given the deterministic trajectory and the model parameters.
    This is accomplished via the model's \code{dmeasure} slot.
    The negative log likelihood is the objective function's value.
  }

  The objective function itself --- in a form suitable for use with \code{\link{optim}}-like optimizers --- is created by a call to \code{traj.match.objfun}.
  Specifically, \code{traj.match.objfun} will return a function that takes a single numeric-vector argument that is assumed to cotain the parameters named in \code{est}, in that order.
}
\value{
  \code{traj.match} returns an object of class \code{traj_matched_pomp}.
  This class inherits from class \code{\link[=pomp-class]{pomp}} and contains the following additional slots:
  \describe{
    \item{transform, est}{
      the values of these arguments on the call to \code{traj.match}.
    }
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

  \code{traj.match.objfun} returns a function suitable for use as an objective function in an \code{\link{optim}}-like optimizer.
}
\examples{
pompExample(ou2)
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

\dontrun{
res1 <- traj.match(
    simdata,
    start=guess.p,
    est=c('alpha.1','alpha.3','alpha.4','x1.0','x2.0','tau'),
    maxit=2000,
    method="nloptr",
    algorithm="NLOPT_LN_SBPLX",
    xtol_rel=1e-5, maxeval=1e5
)

summary(res1)
}

pompExample(ricker)
ofun <- traj.match.objfun(ricker,est=c("r","phi"),transform=TRUE)
optim(fn=ofun,par=c(2,0),method="BFGS")

pompExample(bbs)
## some options are passed to the ODE integrator
ofun <- traj.match.objfun(bbs,est=c("Beta","gamma"),transform=TRUE,hmax=0.001,rtol=1e-6)
optim(fn=ofun,par=c(0,-1),method="Nelder-Mead",control=list(reltol=1e-10))
}
\seealso{
  \code{\link{trajectory}},
  \code{\link{pomp}},
  \code{\link{optim}},
  \code{\link[subplex]{subplex}}
}
\keyword{optimize}
\keyword{ts}
