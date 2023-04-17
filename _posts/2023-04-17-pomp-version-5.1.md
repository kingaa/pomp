---
date: 17 April 2023
layout: pomp
title: pomp version 5.1 released
---
    
**pomp** version 5.1 is now on CRAN and coming soon to [a mirror near you](https://cran.r-project.org/mirrors.html).
This release contains one change that may break existing code.
See below for details.

### Breaking changes

- The `dimnames` attributes of various arrays that appear in **pomp**, including arrays of observables, state variables, parameters, covariates, and so on, have been made uniform.
  In particular, when a dimension of an array corresponds to variables with different names, this dimension is itself named "name".
  Previously, its name was sometimes "variable" and sometimes "parameter".
  This change is meant to streamline interaction with the **tidyverse**.

### New features

- The archiving functions `bake` and `stew` now create the archive directory if it does not already exist.
  A warning is generated when this happens.

### Internal changes

- We no longer import from the superseded package **plyr**.
