---
layout: pomp
title: pomp version 5.4 released
---

**pomp** version 5.4 has been released to CRAN and is making its way to [a mirror near you](https://cran.r-project.org/mirrors.html).

This release contains just one feature enhancement, but it eliminates an annoying warning message on some platforms that inadvertently slipped into the last round of revisions.

### Feature enhancements

The archiving function `stew` now saves timing information by default.

### Bug fixes

A bug resulting in a spurious warning about conversion of `NULL` elements that appears on some platforms has been fixed.
