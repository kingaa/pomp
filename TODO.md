# pomp to-do list

## For pomp 2:

- Put indexing codes into pomp_fun (at C level), and trigger with 'gnsi = TRUE'.
- do we need basic-component arguments at higher levels of the calling hierarchy for documentation purposes?  PROBABLY YES.
- manual pages on the organization of the package
    - estimation algorithms
    - diagnostic tools
    - elementary operations: pfilter, probe, simulate, trajectory
    - workhorses
    - auxiliary functions
    - optimization routines
    - plotting methods
    - parallelization tools
    - reproducibility tools
    - extending the package
- make R-level functions for various distributions and transformations
	- `dmultinom`
	- BetaBinom
	- BetaNegBinom
- remove need to specify paramnames when log,logit,barycentric partrans is given
- `spy` methods for derived objects
- better scheme for indicating derivatives and maps in C snippets
- should the default process model be persistence?
	- what would the corresponding `dprocess` be?
- Kalman filter?
- document and test 'logistic' example

- "guide to upgrading to **pomp** version 2"
- `as.data.frame` and ?? for `listies`
- `coef` returns data frame in some circumstances?
- demonstration of Fisher information via `pfilter` on slice designs?
- new `userdata` function to get void *?
- update 'bbs' example (rename to 'bsflu')?

- perhaps `dinit`?
- change specification of horsemen?
	- new class for specifying precompiled native routines (with its own `PACKAGE` argument)
	- deprecate character specification
	- ~~allow "R snippets": expressions evaluated in proper context?~~
- easier interface for lists of probes in `probe`
- documentation on `mifList`, `pmcmcList`, etc.?
- methods to change data (`obs<-`)?
	- perhaps recognized data variables, states, covariates in calls to `pomp` are replaced?
- put Kalman check in 'gompertz' test: **is it correct?**
- put Kalman check in 'ou2' test
- MCMC proposals as pomp slots?
- probes as pomp slots?
- what does a generic `pomp.fun` interface look like?
- more/better demos and examples
- should parameter transformations allow renaming of variables?
	- would require attention to `rw.sd`, e.g.

## For pomp:

- streamline the R snippets so that ... is unnecessary?
	- probably not worth the trouble, since C snippets are so much faster
- add `include` argument to `pomp`?
- deprecate `obs` and `states` arguments in `simulate`
- graceful stopping for optimizers (at least for `nloptr`)
- `pfilterList` object?
- parameter transformations: put `transform` option into each estimation routine (`spect.match`)
- trap errors for LAPACK
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
