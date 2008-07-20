\name{sobol}
\alias{sobol}
\title{Sobol' low-discrepancy sequence}
\description{
  Generate a data-frame containing a Sobol' low-discrepancy sequence.
}
\usage{
sobol(vars, n)
}
\arguments{
  \item{vars}{Named list of ranges of variables.}
  \item{n}{Number of vectors requested.}
}
\details{}
\value{
  \item{sobol}{Returns a data frame with n `observations' of the
    variables in \code{vars}.}
}
\references{Press et al.}
\author{Aaron A. King (kingaa at umich dot edu)}
\examples{
plot(sobol(vars=list(a=c(0,1),b=c(100,200)),100))
plot(sobol(vars=list(a=c(0,1),b=c(100,200),c=c(-1,1)),100))
}
\keyword{design}
