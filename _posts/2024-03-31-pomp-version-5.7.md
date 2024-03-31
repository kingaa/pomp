---
layout: pomp
title: pomp version 5.7 released
---

**pomp** version 5.7 has been released to CRAN and is making its way to [a mirror near you](https://cran.r-project.org/mirrors.html).

This is a bug-fix release.
Specific changes include:

- With **R** version 4.3.3, changes to the `pbeta` function introduced a bug into **pomp**'s `wquant`.
  This bug has been fixed.
- Changes to some of the low-level C codes in **pomp** resolve issues identified by **rchk** during routine package checking at CRAN.
- Some minor corrections and improvements to the documentation have been made.
