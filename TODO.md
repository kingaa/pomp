-------------

# pomp to-do list

## For pomp:

- Rosenzweig-MacArthur example
- Make certain error checks in inline codes depend on NDEBUG flag
- **R** interface to `eeulermultinom`
- replace `melt` with `data.table::melt`?
- YAML interface
- "Getting Started" vignette:
  - ~~finish discussion of smoothing and hindcasting~~
  - include discussion of forecasting
- ARMA benchmarking (a la **spatPomp**)?
- make `mcap` into a generic, provide a method for data frames and a plotting method

- overriding C snippets should overwrite/delete old C snippets, keep one C snippet file
- `options(pomp_cdir)` should take effect *at run time* and not at initial compile time as currently
- `pomp_fun` should be able to take symbols

- in `pfilter`, perhaps filter trajectory (one sample from the smoothing distribution) should go into `states` slot (??)

- rewrite to use `rmultinom`, which is now part of the **R** API (?)

- more demanding tests of `enkf` and `eakf`
- iterated EnKF
- better scheme for indicating derivatives and maps in C snippets

- support for asymmetric MCMC proposals
- `spy` methods for derived objects
- for log-barycentric transformations, check that parameters are contiguous
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

- demonstration of Fisher information via `pfilter` on slice designs?
- new `userdata` function to get void *?

- easier interface for lists of probes in `probe`
- documentation on `mifList`, `pmcmcList`, etc.?
- methods to change data (`obs<-`)?
	- perhaps recognized data variables, states, covariates in calls to `pomp` are replaced?
- MCMC proposals as pomp slots?
- probes as pomp slots?
- what does a generic `pomp_fun` interface look like?

- should parameter transformations allow renaming of variables?
	- would require attention to `rw_sd`, e.g.

- streamline the R snippets so that ... is unnecessary?
	- probably not worth the trouble, since C snippets are so much faster
- add `include` argument to `pomp`?
- graceful stopping for optimizers (at least for `nloptr`)
- trap errors for LAPACK

- one-point SCQL function for possible use in fitting initial conditions
- partial rejection control for `pfilter`?
- adaptive particle numbers in pfilter (?)

- plugin for adaptive tau leaping algorithm.
- plugin for compartmental models

## Documentation

- FAQ needs update with new issues from @MarieAugerMethe:
	- prediction vs smoothing 
- documentation for existing issues with "todo" label

## Helper packages

- parallel mif-farm and multi-start optimization helpers
- number of processors stored in `bake`, `stew`, `freeze` outputs
- parallel likelihood profiling tool (with MCAP)

## Examples

- SDE examples
- LPA beetle examples
- budmoth examples
- variants on SIR: SEIR, SIRS, ...
- more toy models
- more/better examples
	- need examples of objective function methods
	- `pmcmc` examples

-------------
