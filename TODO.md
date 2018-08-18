# pomp to-do list

- "guide to upgrading to **pomp** version 2"
- all methods need all `pomp` arguments, no?
- make `pomp` return as quickly as possible when no extra work is needed
- all methods use transforms if they exist
- `as.data.frame` and ?? for `listies`
- index by concept using `\concept{}` in the help files
- document user-data functions
- provide full suite of methods for `traj.match.objfun`
- make default initializer into a `pomp_fun` (and obviate the need for the `default_initializer` flag
    - all initializers will require `statenames`!
- matching functions will go away, in favor of objective function constructors:
    - 'probe.match.objfun' looks good: needs testing
    - 'traj.match.objfun' looks good: needs testing
    - 'spect.match.objfun' looks good: needs testing
	- 'nlf' needs to be reworked
- revisit whether `pomp` can handle all `params` chopping and changing
- if `pomp` becomes essentially hidden, where will `zeronames` be documented?
- specific workhorses in individual functions!!
- consolidate examples somehow (it's bad to have so many NULL R files)
- new `userdata` function to get void *?
- update 'bbs' example (rename to 'bsflu')?

## For pomp 2:

- on-the-fly modification of basic components
	- needed for nlf
- better scheme for indicating derivatives and maps in C snippets
-- perhaps `dinit`?
- covariates:
	- ~~covariates provided using a constructor function~~
   	- piecewise constant interpolation of covariates
- change specification of horsemen
	- new class for specifying precompiled native routines (with its own `PACKAGE` argument)
	- deprecate character specification
	- allow "R snippets": expressions evaluated in proper context?
- remove `obs` and `states` arguments in `simulate`? or hide them? `do_simulate` in C?
- easier interface for lists of probes in `probe`
- documentation on `mifList`, `pmcmcList`, etc.
- methods to change data (`obs<-`)?
- refurbish entire test suite
	- perhaps use more, smaller method and function specific tests?
- put Kalman check in 'gompertz' test: **is it correct?**
- put Kalman check in 'ou2' test
- MCMC proposals as pomp slots?
- probes as pomp slots?
- what does a generic `pomp.fun` interface look like?
- more/better demos and examples

## For pomp:

- add `include` argument to `pomp`?
- deprecate `obs` and `states` arguments in `simulate`
- graceful stopping for optimizers (at least for `nloptr`)
- `pfilterList` object?
- parameter transformations: put `transform` option into each estimation routine (`spect.match`)
- trap errors for LAPACK
- better unit tests for `sannbox`
- support for asymmetric MCMC proposals
- one-point SCQL function for possible use in fitting initial conditions
- partial rejection control for `pfilter`?
- adaptive particle numbers in pfilter (?)
- parallel `pfilter` algorithm (?)
- plugin for adaptive tau leaping algorithm.
- plugin for compartmental models

## Documentation

- FAQ needs update with new issues
- Complete rewrite of manual for v. 2

## Helper packages

- parallel mif-farm and multi-start optimization helpers
- number of processors stored in `bake`, `stew`, `freeze` outputs
- parallel likelihood profiling tool

## Examples

- SDE examples
- LPA beetle examples
- budmoth examples

## Tests

### Unreconstructed:

tests/logistic.R
tests/ou2.R
tests/rw2.R

tests/filtfail.R
tests/forecast.R

tests/gillespie2.R
tests/gillespie.R

tests/nlf.R

tests/trajmatch.R

### Groups of related functions:

R/builder.R
R/plugins.R

R/nlf.R
R/nlf_objfun.R
R/nlf_funcs.R
R/nlf_guts.R

R/probe_match.R

R/spect_match.R

R/traj_match.R

### Examples:

- ou2.R
- rw2.R
- gompertz.R
