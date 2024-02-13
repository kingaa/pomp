---
title: Smoothed trajectories via PMCMC
layout: pomp
---

`pmcmc` and `pfilter` now have the capability of saving filtered trajectories.  These can be extracted using the new method `filter.traj`.
The principal use will be in conjunction with `pmcmc`, where, upon convergence to the posterior, samples from the filtered trajectories will be draws from the posterior P[x[1:T] | y[1:T]], i.e., the smoothing distribution.
Thanks to Sebastian Funk for initiating this development!
