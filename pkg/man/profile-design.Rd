\name{profileDesign}
\alias{profileDesign}
\alias{profile.design}
\title{Design matrices for likelihood profile calculations.}
\description{
  \code{profileDesign} generates a data-frame where each row can be used as the starting point for a profile likelihood calculation.
}
\usage{
profileDesign(\dots, lower, upper, nprof)
profile.design(\dots, lower, upper, nprof)
}
\arguments{
  \item{\dots}{
    Specifies the parameters over which to profile.
  }
  \item{lower, upper}{
    Named numeric vectors, specifying the range over which the other parameters are to be sampled.
  }
  \item{nprof}{
    The number of starts per profile point.
  }
}
\value{
  \code{profileDesign} returns a data frame with \code{nprof} points per profile point.
  The other parameters in \code{vars} are sampled using \code{sobol}.
  The second form (\code{profile.design}) is now deprecated.
}
\seealso{\code{\link{sobol}}}
\author{Aaron A. King \email{kingaa at umich dot edu}}
\examples{
## A one-parameter profile design:
x <- profileDesign(p=1:10,lower=c(a=0,b=0),upper=c(a=1,b=5),nprof=20)
dim(x)
plot(x)
## A two-parameter profile design:
x <- profileDesign(p=1:10,q=3:5,lower=c(a=0,b=0),upper=c(b=5,a=1),nprof=20)
dim(x)
plot(x)
}
\keyword{design}
