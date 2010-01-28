\name{LondonYorke}
\alias{LondonYorke}
\docType{data}
\title{Reported cases of chickenpox, measles, and mumps from Baltimore and New York, 1928--1972}
\description{
  \code{LondonYorke} is a data-frame containing the monthly number of reported cases of three childhood diseases from two American cities in the mid-20th century.
}
\usage{data(LondonYorke)}
\examples{
data(LondonYorke)

plot(cases~time,data=LondonYorke,subset=disease=="measles",type='n',main="measles",bty='l')
lines(cases~time,data=LondonYorke,subset=disease=="measles"&town=="Baltimore",col="red")
lines(cases~time,data=LondonYorke,subset=disease=="measles"&town=="New York",col="blue")
legend("topright",legend=c("Baltimore","New York"),lty=1,col=c("red","blue"),bty='n')

plot(cases~time,data=LondonYorke,subset=disease=="chickenpox"&town=="New York",type='l',col="blue",main="chickenpox, New York",bty='l')

plot(cases~time,data=LondonYorke,subset=disease=="mumps"&town=="New York",type='l',col="blue",main="mumps, New York",bty='l')
}
\references{
  W. P. London and J. A. Yorke,
  Recurrent Outbreaks of Measles, Chickenpox and Mumps: I. Seasonal Variation in Contact Rates,
  American Journal of Epidemiology, 98:453--468, 1973.
  }
\seealso{\code{\link{pomp-class}} and the vignettes}
\keyword{datasets}
