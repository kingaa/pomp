\name{Example pomp models}
\alias{pompExample}
\title{Examples of the construction of POMP models.}
\alias{Example pomp models}
\alias{pompExample}
\description{
  \code{pompExample} loads pre-built example \code{pomp} objects.
}
\usage{
pompExample(example, \dots, envir = .GlobalEnv)
}
\arguments{
  \item{example}{
    example to load given as a name or literal character string.
    Evoked without an argument, \code{pompExample} lists all available examples.
  }
  \item{\dots}{
    additional arguments define symbols in the environment within which the example code is executed.
  }
  \item{envir}{
    the environment into which the objects should be loaded.
    If \code{envir=NULL}, then the created objects are returned in a list.
  }
}
\details{
  Directories in the the global option \code{pomp.examples} (set using \code{options()}) are searched for files named \file{example.R}.
  If found, this file will be \code{source}d in a temporary environment.
  Additional arguments to \code{pompExample} define variables within this environment and will therefore be available when the code in \file{example.R} is \code{source}d.
}
\value{
  By default, \code{pompExample} has the side effect of creating one or more \code{pomp} objects in the global workspace.
  If \code{envir=NULL}, there are no side effects; rather, the \code{pomp} objects are returned as a list.
}
\author{Aaron A. King \email{kingaa at umich dot edu}}
\examples{
  pompExample()
  pompExample(euler.sir)
  pompExample("gompertz")
  pompExample(ricker,envir=NULL)
  file.show(system.file("examples/bbs.R",package="pomp"))
}
\seealso{
  \code{\link{blowflies}}, \code{\link{dacca}}, \code{\link{gompertz}},
  \code{\link{ou2}}, \code{\link{ricker}}, \code{\link{rw2}},
  \code{\link{euler.sir}}, \code{\link{gillespie.sir}},
  \code{\link{bbs}}, \code{\link{verhulst}}
}
\keyword{models}
\keyword{datasets}
