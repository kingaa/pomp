---
title: Bug fix for logmeanexp
layout: pomp
---

A bug in `logmeanexp` has been fixed.
Previously, when `se = TRUE`, this function used a delta-method estimate of the variance in `log(mean(exp(x)))` which was accurate when `exp(x)` had small variance but had two problems:
<!--more-->

1. It was overly conservative, being an estimate of the standard error on each element of `x`, as opposed to the standard error on `log(mean(exp(x)))`, as might have been expected
2. It performed poorly when the variance in `exp(x)` was large.
In particular, it returned misleadingly small estimates of the standard error in this case.

The new version uses a jackknife estimate of the variance in `log(mean(exp(x)))`.
Since the jackknife estimate is biased upward, it is still somewhat conservative, but is more robust when the variance of `exp(x)` is large.
