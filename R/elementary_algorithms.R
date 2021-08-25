##' Elementary computations on POMP models.
##'
##' In \pkg{pomp}, elementary algorithms perform POMP model operations.
##' These operations do not themselves estimate parameters, though they may be instrumental in inference methods.
##' 
##' There are six elementary algorithms in \pkg{pomp}:
##' \itemize{
##' \item \code{\link{simulate}} which simulates from the joint distribution of latent and observed variables,
##' \item \code{\link{pfilter}}, which performs a simple particle filter operation,
##' \item \code{\link{wpfilter}}, which performs a weighted particle filter operation,
##' \item \code{\link{probe}}, which computes a suite of user-specified summary statistics to actual and simulated data,
##' \item \code{\link{spect}}, which performs a power-spectral density function computation on actual and simulated data,
##' \item \code{\link{trajectory}}, which iterates or integrates the deterministic skeleton (according to whether the latter is a (discrete-time) map or a (continuous-time) vectorfield.
##' }
##'
##' Help pages detailing each elementary algorithm component are provided.
##'
##' @name elementary algorithms
##' @rdname elementary_algorithms
##' @family elementary algorithms
##' @seealso \link[=basic components]{basic model components},
##' \link[=workhorses]{workhorse functions},
##' \link[=estimation algorithms]{estimation algorithms}.
##' 
NULL
