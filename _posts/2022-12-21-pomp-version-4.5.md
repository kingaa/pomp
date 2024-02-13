---
layout: pomp
title: version 4.5 released to CRAN
---
    
Version 4.5 of **pomp** is now on CRAN and coming soon to [a mirror near you](https://cran.r-project.org/mirrors.html).
This release contains a bug fix and a few changes to the user interface.
Highlights include:

- A bug in trajectory computation for deterministic, discrete-time maps has been fixed.
Thanks to Felicia Magpantay for finding the bug and reporting it!
- One can now extract *weighted* particles from a [`pfilter` computation](https://kingaa.github.io/manuals/pomp/html/pfilter.html).
Before, one could always extract unweighted particles (i.e., particles post resampling).
The impetus for this comes from [Discussion #181](https://github.com/kingaa/pomp/discussions/181).
- Relatedly, the [`saved.states`](https://kingaa.github.io/manuals/pomp/html/saved_states.html) accessor can retrieve the output in a handy data-frame format, and not just in the (somewhat unwieldy) list format as before.
- The [`wquant` function](https://kingaa.github.io/manuals/pomp/html/wquant.html) for estimating quantiles given weighted samples has been changed.
It now uses the Harrell-Davis estimator, which does not seem to be readily available elsewhere, but which should behave well in the primary **pomp** use case---computing quantiles of state distributions from weighted particle ensembles.

In addition, with this release come some small improvements to the [manual](https://kingaa.github.io/manuals/pomp/).
