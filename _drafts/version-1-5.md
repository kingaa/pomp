---
date: 1 June 2016
title: pomp version 1.5 released
layout: pomp
---

### User-visible changes

- A better interface for specifying a model's deterministic skeleton is provided.
  In the `pomp` function, specify `skeleton=map(f,delta.t)` for a discrete-time skeleton (a map) and `skeleton=vectorfield(f)` for a continuous-time skeleton (a vectorfield).
  The old arguments `skeleton.type` and `skelmap.delta.t` are deprecated and will be removed in a future release, at which point the new interface will be mandatory.
- The `method="mif2"` option to `mif` has been removed.
  Use `mif2` instead.
- The `particles` method (rarely if ever used), has been removed to streamline the `mif` codes.
- `mif2` no longer computes filter means of parameters or state variables.
- The minimum version of \R supported is now 3.1.2.

### New features

- The new argument `show` of `pompExamples` allows one to display the example code instead of executing it.
- `init.state` now has the optional argument `nsim`.
  Using this, one can request multiple initial state vectors per parameter vector.

### Documentation improvements

- The `pfilter` help page has been improved.
  Specifically, the discussion of filtering failures is more explicit and easier to find.

### Bug fixes

- A bug associated with the `rw.sd` argument to `mif2` on Windows platforms has been fixed.

### Under the hood

- `pfilter` now uses less memory when it is asked to run at a single point in parameter space.
- The particle filtering algorithms underlying `mif` and `mif2` have now been completely separated from those underlying `pfilter`, a considerable simplification of the codes.

Please see the package [**NEWS**](http://kingaa.github.io/pomp/NEWS.html) for more details.

