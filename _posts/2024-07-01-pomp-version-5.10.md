---
layout: pomp
title: pomp version 5.10 released
---

**pomp** version 5.10 has been released to CRAN and will soon appear on [CRAN mirrors everywhere](https://cran.r-project.org/mirrors.html).
This release contains a change to the `onestep` rprocess simulator and some minor documentation improvements.

### User-visible change

Suppose `P` is a 'pomp' object, with an rprocess component specified as `onestep(f)`, where `f` is a C snippet or **R** function.
Suppose also that `t==time(P,t0=TRUE)`.
In previous package versions, `f` would be executed to go from `t[i]` to `t[i+1]` if and only if `t[i+1] > t[i]`.
As of version 5.9.1, `f` is executed exactly once even if `t[i] == t[i+1]`.
If `f` is written correctly, this change will introduce no error.
However, if the user's code assumes that the time-step `t[i+1]-t[i]` (furnished to `f` as `delta.t` if `f` is an **R** function and `dt` if it is a C snippet) is strictly positive, then this change may introduce errors.

Please contact me via the [package Issues page](https://github.com/kingaa/pomp/issues/) if you notice a change in the behavior of your codes upon update to **pomp** version 5.10.
