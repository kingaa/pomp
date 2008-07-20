\name{mif-methods}
\docType{methods}
\alias{mif-methods}
\alias{logLik,mif-method}
\alias{logLik-mif}
\alias{conv.rec}
\alias{conv.rec,mif-method}
\alias{conv.rec-mif}
\alias{pred.mean}
\alias{pred.mean,mif-method}
\alias{pred.mean-mif}
\alias{pred.var}
\alias{pred.var,mif-method}
\alias{pred.var-mif}
\alias{filter.mean}
\alias{filter.mean,mif-method}
\alias{filter.mean-mif}
\alias{plot-mif}
\alias{plot,mif-method}
\alias{compare.mif}
\title{Methods of the "mif" class}
\description{Methods of the "mif" class.}
\usage{
\S4method{logLik}{mif}(object, \dots)
conv.rec(object, \dots)
\S4method{conv.rec}{mif}(object, pars, \dots)
pred.mean(object, \dots)
\S4method{pred.mean}{mif}(object, pars, \dots)
pred.var(object, \dots)
\S4method{pred.var}{mif}(object, pars, \dots)
filter.mean(object, \dots)
\S4method{filter.mean}{mif}(object, pars, \dots)
\S4method{plot}{mif}(x, y = NULL, \dots)
compare.mif(z)
}
\arguments{
  \item{object}{The \code{mif} object.}
  \item{pars}{Names of parameters.}
  \item{x}{The \code{mif} object.}
  \item{y}{Ignored.}
  \item{z}{A \code{mif} object or list of \code{mif} objects.}
  \item{\dots}{
    Further arguments (either ignored or passed to underlying functions).
  }
}
\section{Methods}{
  \describe{
    \item{conv.rec}{
      \code{conv.rec(object, pars = NULL)} returns the columns of the convergence-record matrix corresponding to the names in \code{pars}.
      By default, all rows are returned.
    }
    \item{logLik}{
      Returns the value in the \code{loglik} slot.
    }
    \item{mif}{
      Re-runs the MIF iterations.
      See the documentation under \cite{Re-running MIF Iterations} below.
    }
    \item{compare.mif}{
      Given a \code{mif} object or a list of \code{mif} objects, \code{compare.mif} produces a set of diagnostic plots.
    }
    \item{plot}{
      Plots a series of diagnostic plots.
      When \code{x} is a \code{mif} object, \code{plot(x)} is equivalent to \code{compare.mif(list(x))}.
    }
    \item{pred.mean}{
      \code{pred.mean(object, pars = NULL)} returns the rows of the prediction-mean matrix corresponding to the names in \code{pars}.
      By default, all rows are returned.
    }
    \item{pred.var}{
      \code{pred.var(object, pars = NULL)} returns the rows of the prediction-variance matrix corresponding to the names in \code{pars}.
      By default, all rows are returned.
    }
    \item{filter.mean}{
      \code{filter.mean(object, pars = NULL)} returns the rows of the filtering-mean matrix corresponding to the names in \code{pars}.
      By default, all rows are returned.
    }
    \item{predvarplot}{
      \code{predvarplot(object, pars = NULL, mean = FALSE, \dots)} produces a plot of the scaled prediction variances for each parameter.
      This can be used to diagnose a good value of the \code{mif} parameters \code{CC} and \code{T0}.
      If used in this way, one should run \code{mif} with \code{Nmif=1} first.
      Additional arguments in \code{\dots} will be passed to the actual plotting function.
    }
    \item{print}{
      Prints a summary of the \code{mif} object.
    }
    \item{show}{
      Displays the \code{mif} object.
    }
    \item{pfilter}{
      See \code{\link{pfilter-mif}}.
    }
  }
}
\references{
  E. L. Ionides, C. Bret{\'o}, & A. A. King,
  Inference for nonlinear dynamical systems,
  Proc. Natl. Acad. Sci. U.S.A., 103:18438--18443, 2006.
}
\author{Aaron A. King (kingaa at umich dot edu)}
\seealso{\code{\link{mif}}, \code{\link{pomp}}, \code{\link{pomp-class}}, \code{\link{pfilter}}}
\keyword{models}
\keyword{ts}
