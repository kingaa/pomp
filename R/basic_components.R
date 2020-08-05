##' Basic POMP model components.
##'
##' Mathematically, the parts of a POMP model include the latent-state process transition distribution, the measurement-process distribution, the initial-state distribution, and possibly a prior parameter distribution.
##' Algorithmically, each of these corresponds to at least two distinct operations.
##' In particular, for each of the above parts, one sometimes needs to make a random draw from the distribution and sometimes to evaluate the density function.
##' Accordingly, for each such component, there are two basic model components, one prefixed by a \sQuote{r}, the other by a \sQuote{d}, following the usual \R convention.
##'
##' In addition to the parts listed above, \pkg{pomp} includes two additional basic model components: the deterministic skeleton, and parameter transformations that can be used to map the parameter space onto a Euclidean space for estimation purposes.
##'
##' There are thus altogether nine \bold{basic model components}:
##' \enumerate{
##' \item \link[=rprocess_spec]{rprocess}, which samples from the latent-state transition distribution,
##' \item \link[=dprocess_spec]{dprocess}, which evaluates the latent-state transition density,
##' \item \link[=rmeasure_spec]{rmeasure}, which samples from the measurement distribution,
##' \item \link[=dmeasure_spec]{dmeasure}, which evaluates the measurement density,
##' \item \link[=prior_spec]{rprior}, which samples from the prior distribution,
##' \item \link[=prior_spec]{dprior}, which evaluates the prior density,
##' \item \link[=rinit_spec]{rinit}, which samples from the initial-state distribution,
##' \item \link[=skeleton_spec]{skeleton}, which evaluates the deterministic skeleton,
##' \item \link[=parameter_trans]{partrans}, which evaluates the forward or inverse parameter transformations.
##' }
##'
##' Each of these can be set or modified in the \code{pomp} constructor function or in any of the \pkg{pomp} \link[=elementary_algorithms]{elementary algorithms} or \link[=estimation_algorithms]{estimation algorithms} using an argument that matches the basic model component.
##' A basic model component can be unset by passing \code{NULL} in the same way.
##'
##' Help pages detailing each basic model component are provided.
##'
##' @name basic_components
##' @rdname basic_components
##' @family implementation_info
##' @seealso \link[=workhorses]{workhorse functions},
##' \link[=elementary_algorithms]{elementary algorithms},
##' \link[=estimation_algorithms]{estimation algorithms}.
##' 
NULL
