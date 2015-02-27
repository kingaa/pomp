\name{Csnippet}
\title{C code snippets for accelerating computations}
\alias{Csnippet}
\alias{Csnippet-class}
\description{
  For including snippets of C code in \code{pomp} objects.
}
\usage{
  Csnippet(text);
}
\arguments{
  \item{text}{
    character; a snippet of C code.
  }
}
\value{
  An object of class \code{Csnippet}.
}
\section{Using C snippets to accelerate computations}{
  From version 0.50, \pkg{pomp} provides a facility whereby users can define their model's components using inline C code.
  The \code{pomp} constructor combines these \code{Csnippet}s into a compilable C file stored in the \R session's temporary directory.
  They are then compiled (via \code{\link[=SHLIB]{R CMD SHLIB}}) into dynamically loadable shared objects.
  This is then dynamically loaded as needed.
}
\author{Aaron A. King \email{kingaa at umich dot edu}}
\seealso{
  \code{\link{pomp}}, \code{\link{plugins}}, and the tutorials on the \href{http://pomp.r-forge.r-project.org}{package website}.
}
\keyword{programming}
\keyword{interface}
