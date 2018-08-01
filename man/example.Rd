\name{Example pomp models}
\title{Examples of the construction of POMP models}
\alias{pompExample}
\alias{Example pomp models}
\alias{pompExample}
\description{
  \code{pompExample} loads pre-built example \sQuote{pomp} objects.
}
\usage{
pompExample(example, \dots, show = FALSE, envir = .GlobalEnv)
}
\arguments{
  \item{example}{
    example to load given as a name or literal character string.
    Evoked without an argument, \code{pompExample} lists all available examples.
  }
  \item{\dots}{
    additional arguments define symbols in the environment within which the example code is executed.
  }
  \item{show}{
    logical; if \code{TRUE}, display, but do not execute, the example \R code.
  }
  \item{envir}{
    the environment into which the objects should be loaded.
    If \code{envir=NULL}, then the created objects are returned in a list.
  }
}
\details{
  Directories listed in the global option \code{pomp.examples} (which can be changed using \code{options()}) are searched for file named \file{<example>.R}.
  If found, this file will be \code{source}d in a temporary environment.
  Additional arguments to \code{pompExample} define variables within this environment and will therefore be available when the code in \file{<example>.R} is \code{source}d.

  The codes that construct these \sQuote{pomp} objects can be found in the \file{examples} directory in the installed package.
  Do \code{system.file("examples",package="pomp"))} to find this directory.
}
\value{
  By default, \code{pompExample} has the side effect of creating one or more objects in the global workspace.
  If \code{envir=NULL}, there are no side effects; rather, the objects are returned as a list.
}
\author{Aaron A. King}
\examples{
  pompExample()
  pompExample(sir)
  pompExample("gompertz")
  pompExample(ricker,envir=NULL)
\dontrun{
  pompExample(bbs,show=TRUE)
}
}
\seealso{
  \code{\link{blowflies}}, \code{\link{dacca}}, \code{\link{gompertz}},
  \code{\link{ou2}}, \code{\link{ricker}}, \code{\link{rw2}},
  \code{\link{sir}}, \code{\link{sir2}}, \code{\link{bbs}}
}
\keyword{models}
\keyword{datasets}
