\name{sannbox}
\alias{sannbox}
\title{Simulated annealing with box constraints.}
\description{
  \code{sannbox} is a straightforward implementation of simulated annealing with box constraints.
}
\usage{
sannbox(par, fn, control = list(), \dots)
}
\arguments{
  \item{par}{
    Initial values for the parameters to be optimized over.
  }
  \item{fn}{
    A function to be minimized, with first argument the vector of parameters over which minimization is to take place.
    It should return a scalar result.
  }
  \item{control}{
    A named list of control parameters.  See \sQuote{Details}.
  }
  \item{\dots}{ignored.}
}
\details{
  The \code{control} argument is a list that can supply any of the following components:
  \describe{
    \item{\code{trace}}{
      Non-negative integer.
      If positive, tracing information on the progress of the optimization is produced.
      Higher values may produce more tracing information.
    }
    \item{\code{fnscale}}{
      An overall scaling to be applied to the value of \code{fn} during optimization.
      If negative, turns the problem into a maximization problem.
      Optimization is performed on \code{fn(par)/fnscale}.
    }
    \item{\code{parscale}}{
      A vector of scaling values for the parameters.
      Optimization is performed on \code{par/parscale} and these should be comparable in the sense that a unit change in any element produces about a unit change in the scaled value.
    }
    \item{\code{maxit}}{
      The total number of function evaluations: there is no other stopping criterion.
      Defaults to \code{10000}.
    }
    \item{\code{temp}}{
      starting temperature for the cooling schedule.
      Defaults to \code{1}.
    }
    \item{\code{tmax}}{
      number of function evaluations at each temperature.
      Defaults to \code{10}.
    }
    \item{\code{candidate.dist}}{
      function to randomly select a new candidate parameter vector.
      This should be a function with three arguments, the first being the current parameter vector, the second the temperature, and the third the parameter scaling.
      By default, \code{candidate.dist} is \preformatted{function(par,temp,scale) rnorm(n=length(par),mean=par,sd=scale*temp)}.
    }
    \item{\code{sched}}{
      cooling schedule.
      A function of a three arguments giving the temperature as a function of iteration number and the control parameters \code{temp} and \code{tmax}.
      By default, \code{sched} is \preformatted{function(k,temp,tmax) temp/log(((k-1)\%/\%tmax)*tmax+exp(1))}.

      Alternatively, one can supply a numeric vector of temperatures.
      This must be of length at least \code{maxit}.
    }
  }
}
\value{
  \code{sannbox} returns a list with components:
  \describe{
    \item{counts}{
      two-element integer vector.
      The first number gives the number of calls made to \code{fn}.
      The second number is provided for compatibility with \code{\link{optim}} and will always be NA.
    }
    \item{convergence}{provided for compatibility with \code{\link{optim}}; will always be 0.}
    \item{final.params}{last tried value of \code{par}.}
    \item{final.value}{value of \code{fn} corresponding to \code{final.params}.}
    \item{par}{best tried value of \code{par}.}
    \item{value}{value of \code{fn} corresponding to \code{par}.}
  }
}
\author{
  Daniel Reuman, Imperial College London and Aaron A. King \email{kingaa at umich dot edu}
}
\seealso{
  \code{\link{traj.match}}, \code{\link{probe.match}}.
}
\keyword{optimize}
