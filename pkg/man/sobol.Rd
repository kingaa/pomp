\name{sobol}
\alias{sobol}
\title{Sobol' low-discrepancy sequence}
\description{
  Generate a data-frame containing a Sobol' low-discrepancy sequence.
}
\usage{
sobol(vars, n)
sobol.design(lower, upper, nseq)
}
\arguments{
  \item{vars}{Named list of ranges of variables.}
  \item{lower,upper}{named numeric vectors giving the lower and upper bounds of the ranges, respectively.}
  \item{n,nseq}{Number of vectors requested.}
}
\value{
  \item{sobol}{
    Returns a data frame with \code{n} \sQuote{observations} of the variables in \code{vars}.
  }
  \item{sobol.design}{
    Returns a data frame with \code{nseq} \sQuote{observations} of the variables over the range specified.
  }
}
\references{
  W. H. Press, S. A. Teukolsky, W. T. Vetterling, \& B. P. Flannery,
  Numerical Recipes in C,
  Cambridge University Press, 1992
}
\seealso{\code{\link{slice.design}}, \code{\link{profile.design}}}
\author{Aaron A. King \email{kingaa at umich dot edu}}
\examples{
plot(sobol(vars=list(a=c(0,1),b=c(100,200)),100))
plot(sobol.design(lower=c(a=0,b=100),upper=c(b=200,a=1),100))
}
\keyword{design}
