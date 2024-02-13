---
layout: pomp
title: Phasing out the particle filter tolerance parameter
---

The new development release, version 2.4.1.3, sets in motion an important change in the behavior of key **pomp** algorithms.
All algorithms based on the particle filter---including `pfilter`, `mif2`, `pmcmc`, and `bsmc2`---are affected.
Each of these algorithms has a tolerance, `tol`, which sets the minimum likelihood distinguishable from zero for the purposes of these methods.
This parameter, the rationale for which has always been practical rather than theoretical, has been present since the earliest days of **pomp**, and there is little evidence for its usefulness.
Worse, as we continue to tackle larger and more complicated problems with **pomp**, we are increasingly encountering situations where the default value of `tol` turns out to be inappropriately high.
For these reasons, we have decided to dispense with it entirely.

Since this change will break some existing code, we will accomplish this in stages.
In a forthcoming release, the new default is set to become zero and ultimately we anticipate a version that entirely removes the option to set a nonzero tolerance.
Accordingly, as of 2.4.1.3, a warning is generated whenever `tol` is nonzero---including for the default value.
To make this annoying warning go away, simply set `tol = 0` in any of the affected algorithms: `pfilter`, `mif2`, `pmcmc`, or `bsmc2`.
To reiterate: the default behavior of these algorithms continues unchanged as of this version, but warnings of forthcoming changes are generated.
These are intended to give users time to adjust their **pomp** usage to `tol = 0` before it becomes, first, the default and, ultimately, mandatory.
