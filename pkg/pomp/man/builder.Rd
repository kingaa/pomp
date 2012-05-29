\name{pompBuilder}
\alias{pompBuilder}
\title{Write, compile, link, and build a pomp object using native codes}
\description{
  \code{pompBuilder} is an EXPERIMENTAL facility for producing compiled \code{pomp} objects.
}
\usage{
pompBuilder(data, times, t0, name, statenames, paramnames, 
            rmeasure, dmeasure, step.fn, step.fn.delta.t,
            skeleton, skeleton.type, skelmap.delta.t = 1,
            parameter.transform, parameter.inv.transform,
            \dots, link = TRUE)
}
\arguments{
  \item{data, times, t0}{
    The data, times, and zero-time.
    See \code{\link{pomp}} for more information.
    \code{data} must be a data-frame.
  }
  \item{name}{
    character; the stem of the name for the files that will be produced.
  }
  \item{statenames, paramnames}{
    names of state-variables and parameters, respectively
  }
  \item{rmeasure, dmeasure}{
    C codes implementing the measurement model
  }
  \item{step.fn, step.fn.delta.t}{
    \code{step.fn} is a C code that implements an Euler step function.
    The Euler time-step is \code{step.fn.delta.t}, which should be a positive number.
  }
  \item{skeleton, skeleton.type, skelmap.delta.t}{
    \code{skeleton} is a C code that implements the deterministic skeleton.
    As in \code{pomp}, \code{skeleton.type} indicates whether the skeleton is a map (discrete-time) or vectorfield (continuous-time).
    If the former, \code{skelmap.delta.t} is the time-step of the map.
  }
  \item{parameter.transform, parameter.inv.transform}{
    optional C codes that implement parameter transformations.
    \code{parameter.transform} maps parameters from the estimation scale to the natural scale;
    \code{parameter.inv.transformation} maps them from the natural scale to the estimation scale.
  }
  \item{\dots}{
    additional arguments are passed to \code{\link{pomp}}
  }
  \item{link}{
    logical; if TRUE, the resulting code will be linked after compilation.
  }
}
\value{
  The constructed \code{pomp} object.
  A side-effect is the writing and compilation of a C code into a shared-object library.
  These files will reside in the current working directory (see \code{\link{getwd}}).
  If \code{pompBuilder} has been called with \code{link=FALSE}, this shared-object library must be linked (see \code{\link{dyn.load}}) before the \code{pomp} object can be used.
}
\details{
  \code{pompBuilder} assumes that files can be written to the current working directory and that shared-object libraries can be compiled and linked, i.e., that \code{R CMD SHLIB} will work.
  This will not typically be the case in out-of-the-box Windows installations.
}
\seealso{
  \code{\link{pomp}} and the demos.
}
\author{Aaron A. King \email{kingaa at umich dot edu}}
\examples{
  ## see the demos 
}
