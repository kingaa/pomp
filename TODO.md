# pomp to-do list

- move all model codes out of the package, except gompertz
- new 'userdata' function to get void *?
- update 'bbs' example (rename to 'bsflu')?

## For pomp 2:

- on-the-fly modification of basic components
- `time` variable could retain its original name?
- `initializer` -> `rinit` (and perhaps `dinit`)
- piecewise constant interpolation of covariates
- methods to change data (`obs<-`)?
- change specification of horsemen
	- new class for specifying precompiled native routines (with its own `PACKAGE` argument)
	- deprecate character specification
	- allow "R snippets": expressions evaluated in proper context?
- remove `obs` and `states` arguments in `simulate`
- remove `mif` and `bsmc` altogether
- better scheme for indicating transformed variables in Csnippet parameter transformations (i.e., use `T` for parameters on transformed scale
- easier interface for lists of probes in `probe`
- better names for SIR examples
- more/better demos and examples
- plugin for adaptive tau leaping algorithm.
- plugin for compartmental models
- refurbish entire test suite

## For pomp:

- add `include` argument to `pomp`?
- all examples and demos should use C snippets (except gompertz?)
- ~~include CWD in search-path for include files?~~
- deprecate `obs` and `states` arguments in `simulate`
- ~~deprecate `mif` and `bsmc`~~
- graceful stopping for optimizers (at least for `nloptr`)
- number of processors stored in `bake`, `stew`, `freeze` outputs
- digest of data stored in `bake`, `stew`, `freeze` outputs
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
