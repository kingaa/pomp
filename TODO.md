# pomp to-do list

## For pomp 2:

- piecewise constant interpolation of covariates
- on-the-fly modification of basic components (pomp 2?)
- `time` variable could retain its original name?
- `initializer` -> `rinit` and perhaps `dinit`
- methods to change data (`obs<-`)
- change specification of horsemen
	- new class for specifying precompiled native routines (with its own `PACKAGE` argument)
	- deprecate character specification
- deprecate `obs` and `states` arguments in `simulate`
- deprecate `mif` and `bsmc`
- better scheme for indicating transformed variables in Csnippet parameter transformations
- easier interface for lists of probes in `probe`
- refurbish entire test suite

## For pomp:

- graceful stopping for optimizers (at least for `nloptr`)
- number of processors stored in `bake`, `stew`, `freeze` outputs
- digest of data stored in `bake`, `stew`, `freeze` outputs
- documentation on `mifList`, `pmcmcList`, etc.
- parameter transformations: put `transform` option into each estimation routine (`spect.match`)
- objective function for spectrum matching
- objective function for NLF
- trap errors for LAPACK
- write Csnippet support for `onestep.dens` and `gillespie.sim` plugins.
- plugin for adaptive tau leaping algorithm.
- plugin for compartmental models
- better unit tests for `sannbox`
- support for asymmetric MCMC proposals
- one-point SCQL function for possible use in fitting initial conditions
- save particle filtering variance?
    Prediction means are optional.
	Only interesting for end-user if one wants to look at residuals.
- partial rejection control for `pfilter`?
- adaptive particle numbers in pfilter (?)
- parallel `pfilter` algorithm (?)

## Examples

- SDE examples
- LPA beetle examples
- budmoth examples
