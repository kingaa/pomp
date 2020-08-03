# pomp to-do list

## For pomp:

- iterated EnKF
- include 'mcap' function
- for log-barycentric transformations, check that parameters are contiguous
- do we need basic-component arguments at higher levels of the calling hierarchy for documentation purposes?  PROBABLY YES.
- more/better examples
	- need examples of objective function methods
	- `pmcmc` examples
- more toy models
- get the metaphors straight! "horses" vs "hitches" vs "workhorses" ("wagons"?)
- manual pages on the organization of the package
    - diagnostic tools
    - auxiliary functions
    - plotting methods
    - parallelization tools
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

- in EnKF and EAKF, allow matrices to depend on parameters

- streamline the R snippets so that ... is unnecessary?
	- probably not worth the trouble, since C snippets are so much faster
- add `include` argument to `pomp`?
- graceful stopping for optimizers (at least for `nloptr`)
- trap errors for LAPACK
- support for asymmetric MCMC proposals
- one-point SCQL function for possible use in fitting initial conditions
- partial rejection control for `pfilter`?
- adaptive particle numbers in pfilter (?)
- plugin for adaptive tau leaping algorithm.
- plugin for compartmental models

## Documentation

- FAQ needs update with new issues from @MarieAugerMethe
- documentation for existing issues with "todo" label

## Helper packages

- parallel mif-farm and multi-start optimization helpers
- parallel `pfilter` algorithm (in **circumstance**)
- number of processors stored in `bake`, `stew`, `freeze` outputs
- parallel likelihood profiling tool

## Examples

- SDE examples
- LPA beetle examples
- budmoth examples
- variants on SIR
