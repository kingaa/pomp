# pomp to-do list

## For pomp 2:

- include 'mcap' function
- add documentation for existing "todo" issues
- for log-barycentric transformations, check that parameters are contiguous
- do we need basic-component arguments at higher levels of the calling hierarchy for documentation purposes?  PROBABLY YES.
- more/better examples
	- need examples of objective function methods
- more toy models
- get the metaphors straight! "horses" vs "hitches" vs "workhorses" ("wagons"?)
- manual pages on the organization of the package
    - estimation algorithms
    - diagnostic tools
    - elementary operations: pfilter, probe, simulate, spect
    - workhorses
    - auxiliary functions
    - plotting methods
    - parallelization tools
    - reproducibility tools
    - extending the package
- make R-level functions for various distributions and transformations
	- BetaBinom
	- BetaNegBinom
- remove need to specify paramnames when log,logit,barycentric partrans is given? HOW?
- `spy` methods for derived objects
- better scheme for indicating derivatives and maps in C snippets
- should the default process model be persistence?
	- what would the corresponding `dprocess` be?
- Kalman filter?

- `melt` for listies
- `coef` returns data frame in some circumstances?
- demonstration of Fisher information via `pfilter` on slice designs?
- new `userdata` function to get void *?
- ~~`simulate` with `format="d"` should return the same thing as `simulate` %>% `as.data.frame`~~

- perhaps `dinit`?
- easier interface for lists of probes in `probe`
- documentation on `mifList`, `pmcmcList`, etc.?
- methods to change data (`obs<-`)?
	- perhaps recognized data variables, states, covariates in calls to `pomp` are replaced?
- put Kalman check in 'ou2' test
- MCMC proposals as pomp slots?
- probes as pomp slots?
- what does a generic `pomp.fun` interface look like?

- should parameter transformations allow renaming of variables?
	- would require attention to `rw.sd`, e.g.

- weighted particle filter
- in EnKF and EAKF, allow matrices to depend on parameters

## For pomp:

- streamline the R snippets so that ... is unnecessary?
	- probably not worth the trouble, since C snippets are so much faster
- add `include` argument to `pomp`?
- ~~deprecate `obs` and `states` arguments in `simulate`~~
- graceful stopping for optimizers (at least for `nloptr`)
- ~~`pfilterList` object~~
- trap errors for LAPACK
- support for asymmetric MCMC proposals
- one-point SCQL function for possible use in fitting initial conditions
- partial rejection control for `pfilter`?
- adaptive particle numbers in pfilter (?)
- parallel `pfilter` algorithm (?)
- plugin for adaptive tau leaping algorithm.
- plugin for compartmental models

## Documentation

- FAQ needs update with new issues from @MarieAugerMethe

## Helper packages

- parallel mif-farm and multi-start optimization helpers
- number of processors stored in `bake`, `stew`, `freeze` outputs
- parallel likelihood profiling tool

## Examples

- SDE examples
- LPA beetle examples
- budmoth examples
