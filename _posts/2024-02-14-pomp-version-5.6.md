---
layout: pomp
title: pomp version 5.6 released
---

**pomp** version 5.6 has been released to CRAN and is making its way to [a mirror near you](https://cran.r-project.org/mirrors.html).

This release contains changes to the C API and some documentation corrections.

- The C-level functions `set_pomp_userdata` and `unset_pomp_userdata` are no longer exported.
- The C-level function `bspline_eval` is no longer exported.
  B-splines are available in C snippets via the functions `bspline_basis_eval`, `bspline_basis_eval_deriv`, `periodic_bspline_basis_eval`, and `periodic_bspline_basis_eval_deriv`.
  See the [documentation of the **pomp** C API](https://kingaa.github.io/pomp/C_API.html#splines) for details.
