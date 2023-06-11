##' Basic POMP model components.
##'
##' Mathematically, the parts of a \acronym{POMP} model include the latent-state process transition distribution, the measurement-process distribution, the initial-state distribution, and possibly a prior parameter distribution.
##' Algorithmically, each of these corresponds to at least two distinct operations.
##' In particular, for each of the above parts, one sometimes needs to make a random draw from the distribution and sometimes to evaluate the density function.
##' Accordingly, for each such component, there are two basic model components, one prefixed by a \sQuote{r}, the other by a \sQuote{d}, following the usual \R convention.
##'
##' In addition to the parts listed above, \pkg{pomp} includes two additional basic model components: the deterministic skeleton, and parameter transformations that can be used to map the parameter space onto a Euclidean space for estimation purposes.
##' There are also basic model components for computing the mean and variance of the measurement process conditional on the latent-state process.
##'
##' There are thus altogether eleven \bold{basic model components}:
##' \enumerate{
##' \item \link[=rprocess specification]{rprocess}, which samples from the latent-state transition distribution,
##' \item \link[=dprocess specification]{dprocess}, which evaluates the latent-state transition density,
##' \item \link[=rmeasure specification]{rmeasure}, which samples from the measurement distribution,
##' \item \link[=emeasure specification]{emeasure}, which computes the conditional expectation of the measurements, given the latent states,
##' \item \link[=vmeasure specification]{vmeasure}, which computes the conditional covariance matrix of the measurements, given the latent states,
##' \item \link[=dmeasure specification]{dmeasure}, which evaluates the measurement density,
##' \item \link[=prior specification]{rprior}, which samples from the prior distribution,
##' \item \link[=prior specification]{dprior}, which evaluates the prior density,
##' \item \link[=rinit specification]{rinit}, which samples from the initial-state distribution,
##' \item \link[=skeleton specification]{skeleton}, which evaluates the deterministic skeleton,
##' \item \link[=parameter_trans]{partrans}, which evaluates the forward or inverse parameter transformations.
##' }
##'
##' Each of these can be set or modified in the \code{pomp} \link[=pomp]{constructor function} or in any of the \pkg{pomp} \link[=elementary algorithms]{elementary algorithms} or \link[=estimation algorithms]{estimation algorithms} using an argument that matches the basic model component.
##' A basic model component can be unset by passing \code{NULL} in the same way.
##'
##' Help pages detailing each basic model component are provided.
##'
##' @name basic components
##' @rdname basic_components
##' @concept basic model components
##' @family implementation information
##' @seealso \link[=workhorses]{workhorse functions},
##' \link[=elementary algorithms]{elementary algorithms},
##' \link[=estimation algorithms]{estimation algorithms}.
##' 
NULL
