---
date: 19 July 2021
layout: pomp
title: new versions of 'bake' and 'stew'
---

The latest development version (3.4.5.0) is available for download on the **pomp** website.
This release includes a major overhaul of the two reproducibility functions, `bake` and `stew`.
These functions are used for caching the results of expensive calculations.
In the past, if the cached calculation was modified, or if any of its dependencies changed, it was necessary to manually delete the cache file.
As of version 3.4.5.0, **pomp** now silently caches information about the calculation and its dependencies.
If these are modified, recomputation is triggered.
See [the package manual](/pomp/manual/bake.html) for more details.

To install from the package repository, do
```
install.packages("pomp",repos="https://kingaa.github.io/")
```
or 
```
devtools::install_github("kingaa/pomp@3.4.5.0")
```
