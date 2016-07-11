### pomp to-do list

- trap errors for LAPACK
- remove deprecated skeleton specification style
- methods to change data (`obs<-`)
- write Csnippet support for `onestep.dens` and `gillespie.sim` plugins.
- plugin for compartmental models
- on-the-fly modification of basic components (pomp 2?)
- SDE examples
- better unit tests for `sannbox`
- easier interface for lists of probes in `probe`
- support for asymmetric MCMC proposals
- documentation on `mifList`, `pmcmcList`, etc.
- `initializer` -> `rinit` and perhaps `dinit`
- objective function for spectrum matching
- objective function for NLF
- one-point SCQL function for possible use in fitting initial conditions
- save particle filtering variance?
    Prediction means are optional.
	Only interesting for end-user if one wants to look at residuals.
- parameter transformations: put `transform` option into each estimation routine (`spect.match`)
- plugin for adaptive tau leaping algorithm.
- LPA beetle examples
- budmoth examples
- partial rejection control for `pfilter`?
- adaptive particle numbers in pfilter (?)
- parallel `pfilter` algorithm (?)
