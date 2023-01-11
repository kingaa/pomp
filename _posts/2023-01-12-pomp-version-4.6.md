---
date: 12 January 2023
layout: pomp
title: version 4.6 released to CRAN
---
    
Version 4.6 of **pomp** is now on CRAN and coming soon to [a mirror near you](https://cran.r-project.org/mirrors.html).
This release contains some changes to the user interface.

### Changes in function names

The biggest change from the user standpoint is that all **pomp** functions with names.that.contain.dots have been deprecated in favor of functions in snake_case.
That is, every function that had a dot (`.`) in its name has been replaced by a function where every dot is replaced by an underscore (`_`).

This is unfortunately necessary to avoid problems with CRAN checks, which (falsely) assume that certain functions with dotted names are S3 methods.
The old function names will continue to work, with a warning.
In a future release, the deprecated functions will be removed.

To help you adapt your code to the new naming convention, you can download and run the [to_snake_case.R](https://kingaa.github.io/scripts/to_snake_case.R) script.
Its usage is straightforward:  
1. Make a directory and copy all files that you wish to edit into it.
2. In an **R** session, source this script.
3. Call the `to_snake_case()` function with the path to your new directory as its sole argument
4. Examine the differences between the files for correctness.
5. Move the new files back into place.

### Other user-visible changes

- The extractor functions `cond_logLik`, `eff_sample_size`, `filter_mean`, `filter_traj`, `forecast`, `pred_mean`, `pred_var`, and `saved_states` now allow you to retrieve their output in a handy data-frame format, and not just in the (somewhat unwieldy) list or array formats as before.
  This is accomplished via the new `format` argument in each of these functions.
  
- The `logmeanexp` function now computes the effective sample size when the option `ess=TRUE` is selected.
  The effective sample size can be useful in determining the reliability of the `logmeanexp` estimate.
  Additionally, `logmeanexp` now returns a fully-named vector when either `se=TRUE` or `ess=TRUE`.

- **pomp** no longer depends on the superseded package **reshape2**.
  Accordingly, the `melt` function---useful for converting arrays and nested lists into data frames---is no longer re-exported from **reshape2**.
  **pomp** now contains a somewhat stricter and more limited version of this useful function.

- The **magrittr** pipe, `%>%`, is no longer re-exported by **pomp**:
  use the native **R** pipe, `|>`, instead.

- The package now requires **R** version 4.1 at least.

### Low-level changes

- The `dimnames` attributes for the arrays computed in `pfilter` and `pmcmc` computations have changed.
  In particular, whereas in previous versions, the `time` dimension was given names that were character strings composed of decimal representations of the time (difficult to work with and prone to roundoff error), the `time` dimension now is not given names.
  That is, the `time` dimension in these arrays can be accessed by position, not by name.
  If you want to match these result to the observation time, use the `format=data.frame` option in the corresponding extractor function.
  
### Bug fixes

- A bug that resulted in the `gompertz` example differing on Windows systems has been fixed.
