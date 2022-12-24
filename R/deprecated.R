##' Deprecated functions
##'
##' These functions are deprecated and will be removed in a future release.
##'
##' @name deprecated
##' @rdname deprecated
##' @keywords internal
##' @include package.R
##' @param ... all arguments are passed to the replacement function.
##'
NULL

##' @rdname deprecated
##' @aliases as.pomp
##' @export
as.pomp <- function (...) {
  .Deprecated("as_pomp")
  as_pomp(...)
}

##' @rdname deprecated
##' @aliases bspline.basis
##' @export
bspline.basis <- function (...) {
  .Deprecated("bspline_basis")
  bspline_basis(...)
}

##' @rdname deprecated
##' @aliases cond.logLik
##' @export
cond.logLik <- function (...) {
  .Deprecated("cond_logLik")
  cond_logLik(...)
}

##' @rdname deprecated
##' @aliases eff.sample.size
##' @export
eff.sample.size <- function (...) {
  .Deprecated("eff_sample_size")
  eff_sample_size(...)
}

##' @rdname deprecated
##' @aliases filter.mean
##' @export
filter.mean <- function (...) {
  .Deprecated("filter_mean")
  filter_mean(...)
}

##' @rdname deprecated
##' @aliases filter.traj
##' @export
filter.traj <- function (...) {
  .Deprecated("filter_traj")
  filter_traj(...)
}

##' @rdname deprecated
##' @aliases mvn.diag.rw
##' @export
mvn.diag.rw <- function (...) {
  .Deprecated("mvn_diag_rw")
  mvn_diag_rw(...)
}

##' @rdname deprecated
##' @aliases mvn.rw
##' @export
mvn.rw <- function (...) {
  .Deprecated("mvn_rw")
  mvn_rw(...)
}

##' @rdname deprecated
##' @aliases mvn.rw.adaptive
##' @export
mvn.rw.adaptive <- function (...) {
  .Deprecated("mvn_rw_adaptive")
  mvn_rw_adaptive(...)
}

##' @rdname deprecated
##' @aliases periodic.bspline.basis
##' @export
periodic.bspline.basis <- function (...) {
  .Deprecated("periodic_bspline_basis")
  periodic_bspline_basis(...)
}

##' @rdname deprecated
##' @aliases pred.mean
##' @export
pred.mean <- function (...) {
  .Deprecated("pred_mean")
  pred_mean(...)
}

##' @rdname deprecated
##' @aliases pred.var
##' @export
pred.var <- function (...) {
  .Deprecated("pred_var")
  pred_var(...)
}

##' @rdname deprecated
##' @aliases probe.acf
##' @export
probe.acf <- function (...) {
  .Deprecated("probe_acf")
  probe_acf(...)
}

##' @rdname deprecated
##' @aliases probe.ccf
##' @export
probe.ccf <- function (...) {
  .Deprecated("probe_ccf")
  probe_ccf(...)
}

##' @rdname deprecated
##' @aliases probe.marginal
##' @export
probe.marginal <- function (...) {
  .Deprecated("probe_marginal")
  probe_marginal(...)
}

##' @rdname deprecated
##' @aliases probe.mean
##' @export
probe.mean <- function (...) {
  .Deprecated("probe_mean")
  probe_mean(...)
}

##' @rdname deprecated
##' @aliases probe.median
##' @export
probe.median <- function (...) {
  .Deprecated("probe_median")
  probe_median(...)
}

##' @rdname deprecated
##' @aliases probe.nlar
##' @export
probe.nlar <- function (...) {
  .Deprecated("probe_nlar")
  probe_nlar(...)
}

##' @rdname deprecated
##' @aliases probe.period
##' @export
probe.period <- function (...) {
  .Deprecated("probe_period")
  probe_period(...)
}

##' @rdname deprecated
##' @aliases probe.quantile
##' @export
probe.quantile <- function (...) {
  .Deprecated("probe_quantile")
  probe_quantile(...)
}

##' @rdname deprecated
##' @aliases probe.sd
##' @export
probe.sd <- function (...) {
  .Deprecated("probe_sd")
  probe_sd(...)
}

##' @rdname deprecated
##' @aliases probe.var
##' @export
probe.var <- function (...) {
  .Deprecated("probe_var")
  probe_var(...)
}

##' @rdname deprecated
##' @aliases rw.sd
##' @export
rw.sd <- function (...) {
  .Deprecated("rw_sd")
  rw_sd(...)
}

##' @rdname deprecated
##' @aliases saved.states
##' @export
saved.states <- function (...) {
  .Deprecated("saved_states")
  saved_states(...)
}

