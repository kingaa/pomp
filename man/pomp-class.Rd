\name{pomp-class}
\docType{class}
\alias{pomp-class}
\alias{plot,pomp-method}
\alias{plot-pomp}
\alias{print,pomp-method}
\alias{print-pomp}
\alias{show,pomp-method}
\alias{show-pomp}
\alias{time,pomp-method}
\alias{time-pomp}
\alias{data.array,pomp-method}
\alias{data.array-pomp}
\alias{data.array}
\alias{coerce,pomp,data.frame-method}

\title{Partially-observed Markov process}
\description{The class \code{pomp} encodes a partially-observed Markov process.}
\section{Objects from the Class}{
  Objects should be created by calls of the function
  \code{pomp}.  See the documentation for \code{\link{pomp}} for
  usage instructions and important warnings.
}
\section{Slots}{
  \describe{
    \item{data}{An array holding the data.  This array is of dimensions
      \code{nobs} x \code{ntimes}, where \code{nobs} is the number of
      observed variables and \code{ntimes} is the number of times at which
      observations were made.}
    \item{times}{The times corresponding to the observations.
      \code{times} must be a strictly increasing numeric vector.}
    \item{t0}{The zero-time.}
    \item{rprocess}{Function of prototype
      \code{rprocess(xstart,times,params,\dots)} which simulates from the
      unobserved process.}
    \item{dprocess}{Function of prototype
      \code{dprocess(x,times,params,log=FALSE,\dots)} which evaluates
      the likelihood of a sequence of consecutive state transitions.}
    \item{rmeasure}{Function of prototype
      \code{rmeasure(x,times,params,\dots)}
      which simulates from the observation process.}
    \item{dmeasure}{Function of prototype
      \code{dmeasure(y,x,params,log=FALSE,\dots)} which gives the
      likelihood of \code{y} given \code{x}.}
    \item{userdata}{A list containing any objects the user desires.  Using
      this mechanism, the user can store additional information necessary
      for the definition of the model.}
  }
}
\section{Methods}{
  \describe{
    \item{plot}{Plots the data.  Additional arguments are passed to the
      low-level plotting routine.}
    \item{print}{Prints the \code{pomp} object in a nice way.}
    \item{show}{Displays the \code{pomp} object.}
    \item{data.array}{\code{data.array(object)} returns the array of
      observations.  \code{data.array(object,vars)} gives just the
      observations of variables \code{vars}.  \code{vars} may specify
      the variables by position or by name.}
    \item{time}{Returns the vector of observation times.}
    \item{coerce}{A \code{pomp} object can be coerced to a data-frame
      via \code{as(object,"data.frame")}.}
    \item{rprocess}{simulates the process model.  See
      \code{\link{rprocess-pomp}}.}
    \item{dprocess}{evaluates the process model density.  See
      \code{\link{dprocess-pomp}}.}
    \item{rmeasure}{simulates the measurement model. See
      \code{\link{rmeasure-pomp}}.}
    \item{dmeasure}{evaluates the measurement-model density.  See
      \code{\link{dmeasure-pomp}}.}
    \item{simulate}{\code{simulate} can be used to simulate state and/or
      observation trajectories.  See documentation under
      \code{\link{simulate-pomp}}.}
  }
}
\author{Aaron A. King (kingaa at umich dot edu)}
\seealso{\code{\link{pomp}},\link{simulate-pomp}}
\keyword{models}
\keyword{ts}
