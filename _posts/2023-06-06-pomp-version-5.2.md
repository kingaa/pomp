---
layout: pomp
title: pomp version 5.2 released
---

**pomp** version 5.2 has been released to CRAN and is making its way to [a mirror near you](https://cran.r-project.org/mirrors.html).
This release contains one bug fix and some internal improvements.
There should be no user-visible changes in **pomp**'s behavior.

Significant changes include:

- A bug in `stew` resulting in namespace conflicts has been repaired.
- **pomp** no longer imports directly from **tidyverse** functions.
  This reduces the list of packages upon which **pomp** depends from 20 to 5.
  There is, however, a new dependency on the **data.table** package.
