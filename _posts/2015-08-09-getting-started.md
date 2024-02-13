---
title: Improved "getting started" tutorial
layout: pomp
---

The [Getting started with `pomp` tutorial](https://kingaa.github.io/pomp/vignettes/getting_started.html) has been substantially improved.
In particular, the tutorial now contains not only an example of the construction of a `pomp` object, but also

- the construction of a likelihood profile using parallelized `mif2` computations,
- the estimation of a posterior probability density using `pmcmc`, including the use of an adaptive proposal distribution (`mvn.rw.adaptive`) and the estimation of smoothed state trajectories (using `filter.traj`).

<!--more-->

Although the tutorial uses as its example a very simple data analysis (one in fact for which far simpler methods than the ones employed would be sufficient), the methods that it illustrates can be readily transposed for dealing with much more challenging problems.
