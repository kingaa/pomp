##' Facilities for making additional information to basic components
##'
##' When \acronym{POMP} basic components need information they can't get from parameters or covariates.
##'
##' It can happen that one desires to pass information to one of the \acronym{POMP} model \dfn{basic components} (see \link[=basic_components]{here for a definition of this term}) outside of the standard routes (i.e., via model parameters or covariates).
##' \pkg{pomp} provides facilities for this purpose.
##' We refer to the objects one wishes to pass in this way as \dfn{user data}.
##'
##' The following will apply to every \link[=basic_components]{basic model component}.
##' For the sake of definiteness, however, we'll use the \code{rmeasure} component as an example.
##' To be even more specific, the measurement model we wish to implement is
##' \preformatted{
##'       y1 ~ Poisson(x1+theta),  y2 ~ Poisson(x2+theta),}
##' where \code{theta} is a parameter.
##' Although it would be very easy (and indeed far preferable) to include \code{theta} among the ordinary parameters (by including it in \code{params}), we will assume here that we have some reason for not wanting to do so.
##'
##' Now, we have the choice of providing \code{rmeasure} in one of three ways:
##' \enumerate{
##'   \item as an \R function,
##'   \item as a C snippet, or
##'   \item as a procedure in an external, dynamically loaded library.
##' }
##' We'll deal with these three cases in turn.
##'
##' @section When the basic component is specified as an \R function:
##' We can implement a simulator for the aforementioned measurement model so: \preformatted{
##'    f <- function (t, x, params, theta, ...) {
##'       y <- rpois(n=2,x[c("x1","x2")]+theta)
##'       setNames(y,c("y1","y2"))
##'    }}
##' So far, so good, but how do we get \code{theta} to this function?
##' We simply provide an additional argument to whichever \pkg{pomp} algorithm we are employing (e.g., \code{\link{simulate}}, \code{\link{pfilter}}, \code{\link{mif2}}, \code{\link{abc}}, etc.).
##' For example:
##' \preformatted{
##'     simulate(..., rmeasure = f, theta = 42, ...)
##' }
##' where the \code{\dots} represent the other \code{simulate} arguments we might want to supply.
##' When we do so, a message will be generated, informing us that \code{theta} is available for use by the \acronym{POMP} basic components.
##' This warning helps forestall accidental triggering of this facility due to typographical error.
##'
##' @section When the basic component is specified via a C snippet:
##' A C snippet implementation of the aforementioned measurement model is:
##' \preformatted{
##'     f <- Csnippet(r"{
##'      double theta = *get_userdata_double("theta");
##'      y1 = rpois(x1+theta); y2 = rpois(x2+theta);
##'     }")}
##' Here, the call to \code{get_userdata_double} retrieves a \emph{pointer} to the stored value of \code{theta}.
##' Note that, by using \R string literals (\code{r"{}"}) we avoid the need to escape the quotes in the C snippet text.
##'
##' It is possible to store and retrieve integer objects also, using \code{get_userdata_int}.
##'
##' One must take care that one stores the user data with the appropriate storage type.
##' For example, it is wise to wrap floating point scalars and vectors with \code{as.double} and integers with \code{as.integer}.
##' In the present example, our call to simulate might look like
##' \preformatted{
##'     simulate(..., rmeasure = f, theta = as.double(42), ...)
##' }
##'
##' Since the two functions \code{get_userdata_double} and \code{get_userdata_int} return pointers, it is trivial to pass vectors of double-precision and integers.
##'
##' A simpler and more elegant approach is afforded by the \code{globals} argument (see below).
##'
##' @section When the basic component is specified via an external library:
##'
##' The rules are essentially the same as for C snippets.
##' \code{typedef} declarations for the \code{get_userdata_double} and \code{get_userdata_int} are given in the \file{pomp.h} header file and these two routines are registered so that they can be retrieved via a call to \code{R_GetCCallable}.
##' See the \href{https://cran.r-project.org/doc/manuals/R-exts.html}{Writing \R extensions manual} for more information.
##'
##' @section Setting \code{globals}:
##'
##' The use of the userdata facilities incurs a run-time cost.
##' It is often more efficient, when using C snippets, to put the needed objects directly into the C snippet library.
##' The \code{globals} argument does this.
##' See the example below.
##'
##' @example examples/userdata.R
##'
##' @name userdata
##' @rdname userdata
##' @include pomp.R pomp_class.R
##' @family implementation information
##'
NULL
