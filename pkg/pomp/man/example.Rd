\name{pompExample}
\alias{pompExample}
\title{Pre-built examples of pomp objects.}
\description{
  \code{pompExample} loads pre-built example \code{pomp} objects.
}
\usage{
pompExample(example, envir = .GlobalEnv)
}
\arguments{
  \item{example}{
    example to load given as a name or literal character string.
    Evoked without an argument, \code{pompExample} lists all available examples.
  }
  \item{envir}{
    the environment into which the objects should be loaded.
  }
}
\value{
  \code{pompExample} has the side effect of creating one or more \code{pomp} objects in the global workspace.
}
\author{Aaron A. King \email{kingaa at umich dot edu}}
\examples{
  pompExample()
  pompExample(euler.sir)
  pompExample("gompertz")
  file.show(system.file("include/pomp.h",package="pomp"))
}
\seealso{
  \code{\link{blowflies}}, \code{\link{dacca}}, \code{\link{gompertz}},
  \code{\link{ou2}}, \code{\link{ricker}}, \code{\link{rw2}},
  \code{\link{euler.sir}}, \code{\link{gillespie.sir}},
  \code{\link{bbs}}, \code{\link{verhulst}}
}
