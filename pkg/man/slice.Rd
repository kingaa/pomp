\name{slice.design}
\alias{slice.design}
\alias{profile.design}
\title{Design matrices for likelihood slices and profiles}
\description{
  \code{slice.design} generates a data-frame representing points taken along one or more slices through a point in a multidimensional space.
  \code{profile.design} generates a data-frame where each row can be used as the starting point for a profile likelihood calculation.
}
\usage{
slice.design(vars, n)
profile.design(\dots, vars, n)
}
\arguments{
  \item{\dots}{
    Specifies the parameters over which to profile.
  }
  \item{vars}{
    For the \code{slice.design} case:
    \code{vars} is a named list of numeric vectors, each of which has length either 1 or 3.
    Variables along which slices are to be taken should have length 3, corresponding to the minimum of the range, central point, and maximum of the range.
    For fixed variables, specify just the value.

    For the \code{profile.design} case:
    Named list of numeric vectors, each of which has length 2, specifying the range over which the parameter is to be sampled.
  }
  \item{n}{
    In the case of \code{slice.design}, the number of points per slice.
    In the case of \code{profile.design}, the number of starts per profile point.
  }
}
\value{
  \code{slice.design} returns a data frame with \code{n} points per slice.
  The column \code{slice} is a factor that tells which slice each point belongs to. 

  \code{profile.design} returns a data frame with \code{n} points per profile point.
  The parameters in \code{vars} are sampled using \code{sobol}.
}
\seealso{\code{\link{sobol}}}
\author{Aaron A. King \email{kingaa at umich dot edu}}
\examples{
## A single 11-point slice through the point c(A=3,B=8,C=0) along the B direction.
x <- slice.design(list(A=3,B=c(0,8,10),C=0),n=11)
dim(x)
plot(x)
## Two slices through the same point along the A and C directions.
x <- slice.design(list(A=c(0,3,5),B=8,C=c(0,0,5)),n=11)
dim(x)
plot(x)

## A one-parameter profile design:
x <- profile.design(p=1:10,vars=list(a=c(0,1),b=c(0,5)),n=20)
dim(x)
plot(x)
## A two-parameter profile design:
x <- profile.design(p=1:10,q=3:5,vars=list(a=c(0,1),b=c(0,5)),n=20)
dim(x)
plot(x)
}
\keyword{design}
