##' The initial-state distribution
##'
##' Specification of rinit
##'
##' To fully specify the unobserved Markov state process, one must give its distribution at the zero-time (\code{t0}).
##' One does this by furnishing a value for the \code{rinit} argument.
##' As usual, this can be provided either as a C snippet or as an \R function.
##' In the former case, bear in mind that:
##' \enumerate{
##'   \item The goal of a this snippet is the construction of a state vector, i.e., the setting of the dynamical states at time \eqn{t_0}{t0}.
##'   \item In addition to the parameters and covariates (if any), the variable \code{t}, containing the zero-time, will be defined in the context in which the snippet is executed.
##'   \item \strong{NB:} The \code{statenames} argument plays a particularly important role when the rinit is specified using a C snippet.
##'    In particular, every state variable must be named in \code{statenames}.
##'    \strong{Failure to follow this rule will result in undefined behavior.}
##'  }
##' \link[=Csnippet]{General rules for writing C snippets can be found here}.
##'
##' If an \R function is to be used, pass
##' \preformatted{
##'    rinit = f
##' }
##' to \code{pomp}, where \code{f} is a function with arguments that can include the initial time \code{t0}, any of the model parameters, and any covariates.
##' As usual, \code{f} may take additional arguments, provided these are passed along with it in the call to \code{pomp}.
##' \code{f} must return a named numeric vector of initial states.
##' It is of course important that the names of the states match the expectations of the other basic components.
##'
##' Note that the state-process \code{rinit} can be either deterministic (as in the default) or stochastic.
##' In the latter case, it samples from the distribution of the state process at the zero-time, \code{t0}.
##'
##' @section Default behavior:
##' By default, \code{pomp} assumes that the initial distribution is concentrated on a single point.
##' In particular, any parameters in \code{params}, the names of which end in \dQuote{\code{_0}} or \dQuote{\code{.0}}, are assumed to be initial values of states.
##' When the state process is initialized, these are simply copied over as initial conditions.
##' The names of the resulting state variables are obtained by dropping the suffix.
##'
##' @name rinit_spec
##' @rdname rinit_spec
##' @family information on model implementation
##'
##' @example examples/rinit_spec.R
##' 
NULL
