\name{pompBuilder}
\alias{pompBuilder}
\title{Write, compile, link, and build a pomp object using native codes}
\description{
  \code{pompBuilder} is an EXPERIMENTAL facility for producing compiled \code{pomp} objects.
}
\usage{
pompBuilder(data, times, t0, name, statenames, paramnames, 
            tcovar, covar,
            rmeasure, dmeasure, step.fn, step.fn.delta.t,
            skeleton, skeleton.type = c("map","vectorfield"),
            skelmap.delta.t = 1,
            parameter.transform, parameter.inv.transform,
            globals, \dots, link = TRUE, save = FALSE)
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
  \item{tcovar, covar}{
    optional.
    \code{covar} is a data-frame containing covariates (variables in columns, timepoints in rows);
    \code{tcovar} is the name of the column containing time.
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
    As in \code{\link{pomp}}, \code{skeleton.type} indicates whether the skeleton is a map (discrete-time) or vectorfield (continuous-time).
    If the former, \code{skelmap.delta.t} is the time-step of the map.
  }
  \item{parameter.transform, parameter.inv.transform}{
    optional C codes that implement parameter transformations.
    \code{parameter.transform} maps parameters from the estimation scale to the natural scale;
    \code{parameter.inv.transformation} maps them from the natural scale to the estimation scale.
  }
  \item{globals}{
    optional.
    C codes here will be included verbatim in the generated code and compiled into the shared-object library.
  }
  \item{\dots}{
    additional arguments are passed to \code{\link{pomp}}
  }
  \item{link}{
    logical; if TRUE, the resulting code will be linked after compilation.
  }
  \item{save}{
    logical; if TRUE, the resulting C code will be saved in the current working directory.
    The shared-object library will also be created in this directory.
    If FALSE, the C code and shared-object library will be created in the temporary directory.
  }
}
\value{
  The constructed \code{pomp} object.
  A side-effect is the writing and compilation of a C code into a dynamics shared object.
  These files will reside in the current working directory (see \code{\link{getwd}}).
  If \code{pompBuilder} has been called with \code{link=FALSE}, this dynamic shared object must be linked (see \code{\link{dyn.load}}) before the \code{pomp} object can be used.

  \strong{Parallel computing caution:}
  note that, since a side-effect of running \code{pompBuilder} is the writing and compilation of a file, caution must be used when running parallel computations.
  Simultaneous execution of \code{pompBuilder} on nodes sharing a common disk will cause all nodes to attempt to write the same file at the same time;
  the results will be undefined.
  If multiple nodes are to work on the same \code{pomp} object built by \code{pompBuilder}, first run \code{pompBuilder} once (perhaps with \code{link=FALSE}) to write the C file and compile it into a dynamic shared object file.
  The \code{pomp} object returned by \code{pompBuilder} can then be shared with all worker nodes.
  Note, however, that all worker nodes must dynamically load (using \code{\link{dyn.load}}) the dynamic shared object file before anything can be done with the \code{pomp} object.
}
\details{
  \code{pompBuilder} with \code{save=TRUE} assumes that files can be written to the current working directory and that dynamic shared objects can be compiled and linked, i.e., that \code{R CMD SHLIB} will work.
  This will not typically be the case in out-of-the-box Windows installations.
}
\seealso{
  \code{\link{pomp}} and the demos.
}
\author{Aaron A. King \email{kingaa at umich dot edu}}
\examples{
  ## see the demos 
}
