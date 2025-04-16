---
layout: pomp
title: pomp version 6.2 released
---

**pomp** version 6.2 has been released to CRAN and will soon appear on [a mirror near you](https://cran.r-project.org/mirrors.html).
This release contains only performance enhancements and internal changes needed to keep up with developments in **R**.

Specifically, the basic Kalman filter (`kalmanFilter`) has been refactored, resulting in a roughly three-fold speed up.

In addition, internal changes have been made in package C codes in anticipation of forthcoming changes to the C API for R.
