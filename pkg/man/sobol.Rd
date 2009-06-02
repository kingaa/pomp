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
\value{
  \item{sobol}{
    Returns a data frame with \code{n} \sQuote{observations} of the variables in \code{vars}.
  }
}
\references{
  W. H. Press, S. A. Teukolsky, W. T. Vetterling, \& B. P. Flannery,
  Numerical Recipes in C,
  Cambridge University Press, 1992
}
\author{Aaron A. King \email{kingaa at umich dot edu}}
\examples{
plot(sobol(vars=list(a=c(0,1),b=c(100,200)),100))
plot(sobol(vars=list(a=c(0,1),b=c(100,200),c=c(-1,1)),100))
}
\keyword{design}
