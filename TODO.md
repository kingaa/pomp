# pomp to-do list

- new `userdata` function to get void *?
- update 'bbs' example (rename to 'bsflu')?
- index by concept using `\concept{}` in the help files

## For pomp 2:

- on-the-fly modification of basic components
	- OK for simulate, probe, pfilter, abc, bsmc2, pmcmc, mif2, kalman, spect
	- needed for probe_match, nlf, spect_match
	- not to be done for trajectory, traj_match
- better scheme for indicating derivatives and maps in C snippets
-- perhaps `dinit`?
- covariates:
	- ~~covariates provided using a constructor function~~
   	- piecewise constant interpolation of covariates
- change specification of horsemen
	- new class for specifying precompiled native routines (with its own `PACKAGE` argument)
	- deprecate character specification
	- allow "R snippets": expressions evaluated in proper context?
- remove `obs` and `states` arguments in `simulate`? or hide? `do_simulate` in C?
- easier interface for lists of probes in `probe`
- methods to change data (`obs<-`)?
- refurbish entire test suite
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
- number of processors stored in `bake`, `stew`, `freeze` outputs
- documentation on `mifList`, `pmcmcList`, etc.
- `pfilterList` object?
- parameter transformations: put `transform` option into each estimation routine (`spect.match`)
- objective function for spectrum matching
- objective function for NLF
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
R/pomp_fun.R

R/probe_match.R

R/minim.R

R/nlf.R
R/nlf_objfun.R
R/nlf_funcs.R
R/nlf_guts.R

R/spect_match.R

R/traj_match.R

### Examples:

- ou2.R
- rw2.R
- gompertz.R
