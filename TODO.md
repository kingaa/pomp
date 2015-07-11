### pomp to-do list

- documentation on 'mifList', 'pmcmcList', etc.
- 'initializer' -> 'rinit' and perhaps 'dinit'
- documentation: "regular parameters" instead of "non-IVP"
- add 'coef' method for 'mifList' objects
- objective function for spectrum matching
- one-point SCQL function for possible use in fitting initial conditions
- save particle filtering variance?
    Prediction means are optional.
	Only interesting for end-user if one wants to look at residuals.
- write Csnippet support for 'onestep.dens' and 'gillespie.sim' plugins.
- parameter transformations: put 'transform' option into each estimation routine (spect.match)
- unit tests for 'sannbox'
- partial rejection control for 'pfilter'?
- add plugin for adaptive tau leaping algorithm.
- add LPA model examples
- SDE examples
- extended Kalman filter
- plugin for compartmental models
- adaptive particle numbers in pfilter (?)
- parallel 'pfilter' algorithm (?)
- 'transform' argument for pmcmc?
    this leads to problems with the prior (?)
