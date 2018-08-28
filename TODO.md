# pomp to-do list

## For pomp 2:

- make R-level functions for various distributions and transformations
- ~~`valid_names` function to check for non-null, non-zero-length names~~
- ~~`compile = FALSE` option for `pomp`~~
- it should be possible to refer to `times` inside `covariate_table`
- do we need basic-component arguments at higher levels of the calling hierarchy for documentation purposes?
- `nsim` argument for `rprior`?
- should parameter transformations allow renaming of variables?
	- would require attention to `rw.sd`, e.g.
- refurbish entire test suite
	- for each example and each demo,
		- tests of each basic component
		- plot of data, simulation, and particle filter (+ probes?)
		- rough likelihood check
	- order of testing development:
		- ~~workhorses~~
		- simulate, trajectory, probe, pfilter
		- abc, pmcmc, mif2, bsmc2, objfun methods
		- test extractor methods, etc. as appropriate
		- helper functions get their own test scripts
- ~~we don't need `ncovars` in workhorse C function calls~~
- combine documentation: (dmeasure,rmeasure), (dprocess,rprocess)?
- regularize use of `pStop` and `tryCatch` on *.internals
- index by concept using `\concept{}` in the help files
	- extending \pkg{pomp}
	- low-level interface (incl. `hitch` and workhorses)
	- 
- better scheme for indicating derivatives and maps in C snippets
- should the default process model be persistence?
	- what would the corresponding `dprocess` be?
- ~~revisit use of `enum`~~
- covariates:
	- ~~covariates provided using a constructor function~~
   	- piecewise constant interpolation of covariates
- ~~document default workhorses~~
- ~~get rid of 'times' as dimnames~~
- ~~`do_simulate` should handle parameter column names~~
- ~~`bsmc2` should only take vector `params` and rely on `rprior`~~
- ~~`spect` should call `do_simulate`~~
- ~~change `as.data.frame` argument to `format` in `trajectory`~~
- provide full suite of methods for `traj.match.objfun`
- ~~revisit default initializer~~
- ~~should unspecified basic components return NAs?~~

- "guide to upgrading to **pomp** version 2"
- ~~take away option of passing `params` as a matrix to `pfilter`~~
- ~~all methods use transforms if they exist~~
- `as.data.frame` and ?? for `listies`
- ~~document user-data functions~~
- matching functions will go away, in favor of objective function constructors:
    - 'probe.match.objfun' looks good: needs testing
    - 'traj.match.objfun' looks good: needs testing
    - 'spect.match.objfun' looks good: needs testing
	- 'nlf' needs to be reworked
- ~~revisit whether `pomp` can handle all `params` chopping and changing~~
- if `pomp` becomes essentially hidden, where will `zeronames` be documented?
- ~~put all NLF codes into one file~~
- import and re-export `reshape2::melt`?
- demonstration of Fisher information via `pfilter` on slice designs?
- new `userdata` function to get void *?
- update 'bbs' example (rename to 'bsflu')?

- on-the-fly modification of basic components
	- needed for nlf
- perhaps `dinit`?
- change specification of horsemen?
	- new class for specifying precompiled native routines (with its own `PACKAGE` argument)
	- deprecate character specification
	- allow "R snippets": expressions evaluated in proper context?
- ~~remove `obs` and `states` arguments in `simulate`? or hide them?~~
- ~~`do_simulate` in C~~
- easier interface for lists of probes in `probe`
- documentation on `mifList`, `pmcmcList`, etc.?
- methods to change data (`obs<-`)?
	- perhaps recognized data variables, states, covariates in calls to `pomp` are replaceed?
- put Kalman check in 'gompertz' test: **is it correct?**
- put Kalman check in 'ou2' test
- MCMC proposals as pomp slots?
- probes as pomp slots?
- what does a generic `pomp.fun` interface look like?
- more/better demos and examples
- document Beta-binomial and other distributions

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
