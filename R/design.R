##' Design matrices for pomp calculations
##'
##' These functions are useful for generating designs for the exploration of parameter space.
##'
##' @name design
##' @rdname design
##'
##' @author Aaron A. King
##'
##' @references
##' P. Bratley and B. L. Fox, Algorithm 659 Implementing Sobol's quasirandom
##' sequence generator, ACM Trans. Math. Soft. 14, 88--100, 1988.
##'
##' S. Joe and F. Y. Kuo, Remark on algorithm 659: Implementing Sobol's
##' quasirandom sequence generator ACM Trans. Math. Soft 29, 49--57, 2003.
##'
##' Steven G. Johnson, The \pkg{NLopt} nonlinear-optimization package,
##' \url{http://ab-initio.mit.edu/nlopt}
##'
##' @keywords design
##' @example examples/design.R
##' @param ...
##' In \code{profileDesign}, additional arguments specify the parameters over which to profile and the values of these parameters.
##' In \code{sliceDesign}, additional numeric vector arguments specify the locations of points along the slices.
##' 
NULL
