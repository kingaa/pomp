##' B-spline bases
##'
##' These functions generate B-spline basis functions.
##' \code{bspline_basis} gives a basis of spline functions.
##' \code{periodic_bspline_basis} gives a
##' basis of periodic spline functions.
##'
##' @name bsplines
##' @rdname bsplines
##' @concept splines
##' @family interpolation
##'
##' @param x Vector at which the spline functions are to be evaluated.
##' @param rg numeric of length 2; range of the B-spline basis.
##' To be properly specified, we must have \code{rg[1] < rg[2]}.
##' @param nbasis The number of basis functions to return.
##' @param degree Degree of requested B-splines.
##' @param period The period of the requested periodic B-splines.
##' @param deriv The order of the derivative required.
##' @param names optional; the names to be given to the basis functions.  These
##' will be the column-names of the matrix returned.  If the names are
##' specified as a format string (e.g., "basis\%d"), \code{\link{sprintf}} will
##' be used to generate the names from the column number.  If a single
##' non-format string is specified, the names will be generated by
##' \code{\link{paste}}-ing \code{name} to the column number.  One can also
##' specify each column name explicitly by giving a length-\code{nbasis} string
##' vector.  By default, no column-names are given.
##'
##' @return
##' \item{bspline_basis}{ Returns a matrix with \code{length(x)} rows
##' and \code{nbasis} columns.  Each column contains the values one of the
##' spline basis functions.}
##' \item{periodic_bspline_basis}{ Returns a matrix with \code{length(x)} rows
##' and \code{nbasis} columns.  The basis functions returned are periodic with
##' period \code{period}.}
##' If \code{deriv>0}, the derivative of that order of each of the corresponding spline basis functions are returned.
##'
##' @section C API:
##' Access to the underlying C routines is available: see
##' \href{https://kingaa.github.io/pomp/C_API.html}{the \pkg{pomp} C API document}.
##' for definition and documentation of the C API.
##'
##' @author Aaron A. King
##'
##' @keywords smooth
##' @examples
##'
##' x <- seq(0,2,by=0.01)
##' y <- bspline_basis(x,degree=3,nbasis=9,names="basis")
##' matplot(x,y,type='l',ylim=c(0,1.1))
##' lines(x,apply(y,1,sum),lwd=2)
##'
##' x <- seq(-1,2,by=0.01)
##' y <- periodic_bspline_basis(x,nbasis=5,names="spline%d")
##' matplot(x,y,type='l')
##'
NULL

##' @rdname bsplines
##' @export
bspline_basis <- function (x, nbasis, degree = 3, deriv = 0, names = NULL, rg = range(x)) {
  rg <- as.numeric(rg)
  y <- tryCatch(
    .Call(P_bspline_basis,rg,x,nbasis,degree,deriv),
    error = function (e) {
      pStop(who="bspline_basis",conditionMessage(e))
    }
  )
  if (deriv > degree)
    pWarn("returning 0 since ",sQuote("deriv")," > ",sQuote("degree"))
  if (!is.null(names)) {
    if (length(names)==1) {
      if (!grepl("%",names)) {
        names <- paste0(names,".%d")
      }
      colnames(y) <- sprintf(names,seq_len(nbasis))
    } else if (length(names)==nbasis) {
      colnames(y) <- names
    } else {
      pStop(sQuote("names")," must be of length 1 or ",nbasis)
    }
  }
  y
}

##' @rdname bsplines
##' @export
periodic_bspline_basis <- function (x, nbasis, degree = 3, period = 1,
  deriv = 0, names = NULL) {
  y <- tryCatch(
    .Call(P_periodic_bspline_basis,x,nbasis,degree,period,deriv),
    error = function (e) {
      pStop(who="periodic_bspline_basis",conditionMessage(e))
    }
  )
  if (deriv > degree)
    pWarn("returning 0 since ",sQuote("deriv")," > ",sQuote("degree"))
  if (!is.null(names)) {
    if (length(names)==1) {
      if (!grepl("%",names)) {
        names <- paste0(names,".%d")
      }
      colnames(y) <- sprintf(names,seq_len(nbasis))
    } else if (length(names)==nbasis) {
      colnames(y) <- names
    } else {
      pStop(sQuote("names")," must be of length 1 or ",nbasis)
    }
  }
  y
}
