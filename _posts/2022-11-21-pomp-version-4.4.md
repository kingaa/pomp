---
date: 21 November 2022
layout: pomp
title: version 4.4 released to CRAN
---
    
Version 4.4 of **pomp** is now on CRAN and coming soon to [a mirror near you](https://cran.r-project.org/mirrors.html).
This release contains a few changes to the user interface.
Highlights include:

- The new function `wquant` computes weighted quantiles.
- The precise order in which the `pmcmc` function computes prior and likelihood of furnished and proposed parameters has changed slightly.
  This prevents proposals that are incompatible with the prior from being passed to `pfilter` and forestalls an associated class of errors. 
  As a consequence of this change, `pmcmc` computations using the new version will differ very slightly from previous computations, even with the RNG seed fixed at its previous values.
- There is now a C interface to the `bspline_eval` function.
  See the **pomp** C API documentation for details.
- The `bspline.basis` function now takes the optional argument `rg` which allows one to specify the range over which the basis will be constructed.
  By default, this is `range(x)`, which agrees with the behavior in earlier versions.

In addition, there are some small improvements to the [manual](https://kingaa.github.io/manuals/pomp/).
