\name{sliceDesign}
\alias{sliceDesign}
\alias{slice.design}
\title{Design matrices for likelihood slices.}
\description{
  \code{sliceDesign} generates a data-frame representing points taken along one or more slices through a point in a multidimensional space.
}
\usage{
sliceDesign(center, ...)
slice.design(center, ...)
}
\arguments{
  \item{center}{
    \code{center} is a named numeric vectors specifying the point through which the slice(s) is (are) to be taken.
  }
  \item{\dots}{
    Additional numeric vector arguments specify the slices.
  }
}
\value{
  \code{sliceDesign} returns a data frame with one row per point along a slice.
  The column \code{slice} is a factor that tells which slice each point belongs to. 
  The equivalent form (\code{slice.design}) is now deprecated.
}
\seealso{\code{\link{profile.design}}}
\author{Aaron A. King \email{kingaa at umich dot edu}}
\examples{
## A single 11-point slice through the point c(A=3,B=8,C=0) along the B direction.
x <- sliceDesign(center=c(A=3,B=8,C=0),B=seq(0,10,by=1))
dim(x)
plot(x)
## Two slices through the same point along the A and C directions.
x <- sliceDesign(c(A=3,B=8,C=0),A=seq(0,5,by=1),C=seq(0,5,length=11))
dim(x)
plot(x)
}
\keyword{design}
