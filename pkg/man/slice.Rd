\name{slice.design}
\alias{slice.design}
\title{Slices through multidimensional parameter space}
\description{
  Generate a data-frame representing points taken along one or more slices through a point in a multidimensional space.
}
\usage{
slice.design(vars, n)
}
\arguments{
  \item{vars}{
    Named list of numeric vectors, each of which has length either 1 or 3.
    Variables along which slices are to be taken should have length 3, corresponding to the minimum of the range, central point, and maximum of the range.
    For fixed variables, specify just the value.
  }
  \item{n}{Number of points per slice.}
}
\value{
  \code{slice.design} returns a data frame with \code{n} points per slice.
  The column \code{slice} is a factor that tells which slice each point belongs to. 
}
\author{Aaron A. King \email{kingaa at umich dot edu}}
\examples{
## A single 11-point slice through the point c(A=3,B=8,C=0) along the B direction.
slice.design(list(A=3,B=c(0,8,10),C=0),n=11)
## Two slices through the same point along the A and C directions.
slice.design(list(A=c(0,3,5),B=8,C=c(0,0,5)),n=11)
}
\keyword{design}
