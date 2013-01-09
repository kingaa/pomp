\name{pompExample}
\alias{pompExample}
\title{Pre-built examples of pomp objects.}
\description{
  \code{pompExample} loads pre-built example \code{pomp} objects.
}
\usage{
pompExample(example)
}
\arguments{
  \item{example}{
    character; name of the example to load.
    Evoked without an argument, \code{pompExample} lists all available examples.
  }
}
\value{
  \code{pompExample} has the side effect of creating one or more \code{pomp} objects in the global workspace.
}
\author{Aaron A. King \email{kingaa at umich dot edu}}
\examples{
  pompExample()
  pompExample("euler.sir")
}
