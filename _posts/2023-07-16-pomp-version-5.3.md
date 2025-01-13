---
layout: pomp
title: pomp version 5.3 released
---

**pomp** version 5.3 has been released to CRAN and is making its way to [a mirror near you](https://cran.r-project.org/mirrors.html).

This release contains a number of feature enhancements and a few bug fixes.

### Feature enhancements

- This version adds the new basic component `dinit` for evaluating the probability density function of the initial-state distribution.  There is a corresponding workhorse function, `dinit()`.
	
- The constructor `pomp` now takes the optional argument `nstatevars`, which can be used to increase the dimension of the latent state-vectors created by `rinit`.  By default, `nstatevars = length(statenames)`, and `nstatevars` can only be used to _increase_, not to decrease, the dimension of the latent state process.  Moreover, `nstatevars` has no effect if the `rinit` basic component is furnished as an R function.

- The data frames returned by `cond_logLik` and `eff_sample_size` when `format="data.frame"` have been improved.  In particular, they contain (as variable `time`) the times (rather than the index of the time vector as before).
	
- The new option `on_load` allows one to specify a C snippet that will be executed at the time the C snippet library is loaded.

- This release includes some cosmetic changes to the report generated by `spy`.
  
- The manual pages have been reorganized, with improved cross-linking.

- `reulermultinom` now returns `NA` rather than `NaN`, in keeping with the behavior of `rbinom`.  Thanks to John Drake for calling attention to this issue.

### Internal changes
	
- Internally, the `userdata` are made available via calls to `pompLoad` rather than within each call to a `pomp` workhorse.
	
- Name checking on internally-created vectors and arrays is less strict than previously.  In particular, it is possible to have variables without names (i.e., `""`).

### Bug fixes

- A bug in `filter_traj` etc. arising when not all state variables have names was fixed.