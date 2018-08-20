##' Historical childhood disease incidence data
##'
##' \code{LondonYorke} is a data frame containing the monthly number of
##' reported cases of chickenpox, measles, and mumps from two American cities
##' (Baltimore and New York) in the mid-20th century (1928--1972).
##'
##' \code{ewmeas} and \code{ewcitmeas} are data frames containing weekly
##' reported cases of measles in England and Wales.  \code{ewmeas} records the
##' total measles reports for the whole country, 1948--1966.  One questionable
##' data point has been replaced with an NA.  \code{ewcitmeas} records the
##' incidence in seven English cities 1948--1987.  These data were kindly
##' provided by Ben Bolker, who writes: \dQuote{Most of these data have been
##' manually entered from published records by various people, and are prone to
##' errors at several levels. All data are provided as is; use at your own
##' risk.}
##'
##'
##' @name measles
##' @rdname measles
##' @aliases LondonYorke ewmeas ewcitmeas
##' @docType data
##' @references W. P. London and J. A. Yorke, Recurrent Outbreaks of Measles,
##' Chickenpox and Mumps: I. Seasonal Variation in Contact Rates, American
##' Journal of Epidemiology, 98:453--468, 1973.
##' @keywords datasets
##' @family pomp examples
##' @examples
##'
##' plot(cases~time,data=LondonYorke,subset=disease=="measles",type='n',main="measles",bty='l')
##' lines(cases~time,data=LondonYorke,subset=disease=="measles"&town=="Baltimore",col="red")
##' lines(cases~time,data=LondonYorke,subset=disease=="measles"&town=="New York",col="blue")
##' legend("topright",legend=c("Baltimore","New York"),lty=1,col=c("red","blue"),bty='n')
##'
##' plot(
##'      cases~time,
##'      data=LondonYorke,
##'      subset=disease=="chickenpox"&town=="New York",
##'      type='l',col="blue",main="chickenpox, New York",
##'      bty='l'
##'     )
##'
##' plot(
##'      cases~time,
##'      data=LondonYorke,
##'      subset=disease=="mumps"&town=="New York",
##'      type='l',col="blue",main="mumps, New York",
##'      bty='l'
##'     )
##'
##' plot(reports~time,data=ewmeas,type='l')
##'
##' plot(reports~date,data=ewcitmeas,subset=city=="Liverpool",type='l')
##'
NULL
