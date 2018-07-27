# pomp to-do list

- new `userdata` function to get void *?
- update 'bbs' example (rename to 'bsflu')?
- ~~make new header file for package developers~~
- ~~move all model codes out of the package, except gompertz and ou2~~
- ~~export `hitch` function for compiling `pomp.fun`s~~
- index by concept using `\concept{}` in the help files

## For pomp 2:

- on-the-fly modification of basic components
- `time` variable could retain its original name?
- `initializer` -> `rinit` (and perhaps `dinit`)
- `conv.rec` -> `traces` (with deprecation warning)
- MCMC proposals as pomp slots?
- piecewise constant interpolation of covariates
- methods to change data (`obs<-`)?
- remove `$` methods and replace with proper extractor functions
- change specification of horsemen
	- new class for specifying precompiled native routines (with its own `PACKAGE` argument)
	- deprecate character specification
	- allow "R snippets": expressions evaluated in proper context?
- remove `obs` and `states` arguments in `simulate`? or hide? `do_simulate` in C?
- ~~remove `mif` and `bsmc` altogether~~
- better scheme for indicating transformed variables in C snippet parameter transformations (i.e., use `T` for parameters on transformed scale
- better scheme for indicating derivatives and maps in C snippets
- easier interface for lists of probes in `probe`
- better names for SIR examples
- more/better demos and examples
- plugin for adaptive tau leaping algorithm.
- plugin for compartmental models
- refurbish entire test suite

## For pomp:

- add `include` argument to `pomp`?
- ~~include CWD in search-path for include files?~~
- deprecate `obs` and `states` arguments in `simulate`
- ~~deprecate `mif` and `bsmc`~~
- graceful stopping for optimizers (at least for `nloptr`)
- number of processors stored in `bake`, `stew`, `freeze` outputs
- documentation on `mifList`, `pmcmcList`, etc.
- `pfilterList` object?
- parameter transformations: put `transform` option into each estimation routine (`spect.match`)
- objective function for spectrum matching
- objective function for NLF
- trap errors for LAPACK
- ~~write Csnippet support for `onestep.dens` and `gillespie.sim` plugins.~~
- better unit tests for `sannbox`
- support for asymmetric MCMC proposals
- one-point SCQL function for possible use in fitting initial conditions
- save particle filtering variance?
  Prediction means are optional.
	Only interesting for end-user if one wants to look at residuals.
- partial rejection control for `pfilter`?
- adaptive particle numbers in pfilter (?)
- parallel `pfilter` algorithm (?)

## Helper packages

- parallel mif-farm and multi-start optimization helpers
- parallel likelihood profiling tool

## Examples

- SDE examples
- LPA beetle examples
- budmoth examples

## Tests

### Unreconstructed:

tests/bsplines.R
tests/design.R

tests/prior.R
tests/procmeas.R
tests/proposals.R

tests/logistic.R
tests/ou2.R
tests/rw2.R
tests/sir2.R

tests/dimchecks.R
tests/dp.R
tests/fhn.R
tests/filtfail.R
tests/forecast.R

tests/getting_started.R

tests/gillespie2.R
tests/gillespie.R
tests/issue56.R

tests/kalman.R

tests/mif2-index.R
tests/mif2.R

tests/nlf.R

tests/pomp.R

tests/probe.R

tests/sannbox.R

tests/spect.R

tests/steps.R

tests/synlik.R

tests/trajmatch.R

### Groups of related functions:

R/generics.R

R/parmat.R
R/sobol.R
R/profile_design.R
R/slice_design.R
R/bsplines.R
R/logmeanexp.R

R/builder.R
R/load.R
R/plugins.R
R/pomp_class.R
R/pomp_fun.R
R/pomp_methods.R
R/pomp.R

R/initstate_pomp.R
R/dprior_pomp.R
R/dprocess_pomp.R
R/rprior_pomp.R
R/rprocess_pomp.R
R/trajectory_pomp.R

R/basic_probes.R
R/probe_match.R
R/probe.R

R/kalman_methods.R
R/kalman.R

R/mif2_methods.R
R/mif2.R

R/minim.R

R/nlf_funcs.R
R/nlf_guts.R
R/nlf_objfun.R
R/nlf.R

R/plot_pomp.R

R/proposals.R

R/covmat.R
R/sannbox.R

R/simulate_pomp.R

R/spect_match.R
R/spect.R

R/traj_match.R

### Examples:

- ou2.R
- rw2.R

### demos:

- gompertz.R
- logistic.R
- rw2.R
- sir.R
