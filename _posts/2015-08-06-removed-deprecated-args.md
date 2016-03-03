---
date: 6 Aug 2015
title: Deprecated arguments removed
layout: pomp
---

A number of minor changes to **pomp** have been made since the last blog post.
<!--more-->
Specifically, in **pomp** version 0.77-1:

- The deprecated `transform.params` argument to `nlf` has been removed.
- The deprecated `pars` and `rw.sd` arguments to `abc` and `pmcmc` have been removed.
- The deprecated `data.array` method has been removed.
  Use `obs` instead.
- The default value of `max.fail` for `pmcmc` is now `Inf`.
  That is, by default `pmcmc` will not halt (but will issue a warning) when filtering failure occurs (i.e., no particles are consistent with the data, as defined by the likelihood threshold `tol`).
