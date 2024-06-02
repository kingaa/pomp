---
layout: pomp
title: pomp version 5.9 released
---

**pomp** version 5.9 has been released to CRAN and will soon appear on [a mirror near you](https://cran.r-project.org/mirrors.html).
This release contains an important change to the user interface, as well as a bug fix and one new experimental feature.

## User-visible change

As documented in the manual, [basic model components](https://kingaa.github.io/manuals/pomp/html/basic_components.html) have access to parameters, covariates, state variables, and in some cases, the data.
To make additional elements available to the basic model components, **pomp** provides the ["userdata" facility](https://kingaa.github.io/manuals/pomp/html/userdata.html).
From this version, such elements should be furnished to **pomp** functions in a list via the argument `userdata`.
For the present, the old behavior will still work, but will generate a warning.
In a future release, this will become an error.

## Bug fix
	
A bug in covariate-table extrapolation for the case `order="constant"` has been fixed.

## Feature enhancement

A new experimental function `repair_lookup_table` is provided to help eliminate unnecessary warnings about extrapolation.
One can use it "repair" a covariate table that does not span the full temporal range of observations.
This is simply a matter of performing the extrapolations and adding them to the lookup table.
See the [manual](https://kingaa.github.io/manuals/pomp/html/covariate_table.html) for more information and an example.

## Other changes

Functions with names that include dots (`.`), which have been defunct since version 4.7, have now been expunged entirely.
