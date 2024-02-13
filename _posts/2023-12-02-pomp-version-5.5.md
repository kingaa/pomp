---
layout: pomp
title: pomp version 5.5 released
---

**pomp** version 5.5 has been released to CRAN and is making its way to [a mirror near you](https://cran.r-project.org/mirrors.html).

This release contains a few feature enhancements and some changes to the internals.

### Feature enhancements

- There is now access to ordinary B-spline basis functions (and their derivatives) at the C snippet level.
  It will still almost always be preferable to construct a spline basis once and pass it to **pomp** functions using `covar=covariate_table(...)`, but this functionality may occasionally be useful.
- `bspline_eval` is now deprecated as part of the **pomp** C API.
- Stateful objective functions (class `objfun`) now behave as if they are objects of class `pomp`.
  This makes it easier to work with them in plotting results, extracting information from the fitted model object, etc.

### Internal changes

- The C-level functions `set_pomp_userdata` and `unset_pomp_userdata` have been deprecated and will be removed in a future release.
  Using the latter now generates a warning that can be safely ignored.
