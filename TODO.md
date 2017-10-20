### pomp to-do list

- on-the-fly modification of basic components (pomp 2?)
- `time` variable could retain its original name?
- `initializer` -> `rinit` and perhaps `dinit`
- come up with better scheme for indicating transformed variables in Csnippet parameter transformations
- parameter transformations: put `transform` option into each estimation routine (`spect.match`)
- change specification of horsemen
	- new class for specifying precompiled native routines (with its own `PACKAGE` argument)
	- deprecate character specification
- objective function for spectrum matching
- objective function for NLF
- trap errors for LAPACK
- methods to change data (`obs<-`)
- write Csnippet support for `onestep.dens` and `gillespie.sim` plugins.
- plugin for compartmental models
- SDE examples
- better unit tests for `sannbox`
- easier interface for lists of probes in `probe`
- support for asymmetric MCMC proposals
- documentation on `mifList`, `pmcmcList`, etc.
- one-point SCQL function for possible use in fitting initial conditions
- save particle filtering variance?
    Prediction means are optional.
	Only interesting for end-user if one wants to look at residuals.
- plugin for adaptive tau leaping algorithm.
- LPA beetle examples
- budmoth examples
- partial rejection control for `pfilter`?
- adaptive particle numbers in pfilter (?)
- parallel `pfilter` algorithm (?)
