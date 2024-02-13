---
layout: pomp
title: new document defining the pomp C API
---

**pomp** provides C entry points to a number of facilities for model specification.
The new [C API document](https://kingaa.github.io/pomp/C_API.html) describes these.

In particular, **pomp** provides C-level access to several different probability distributions of use in modeling, facilities for working with splines, and parameter transformations.

The [final section of the document](https://kingaa.github.io/pomp/C_API.html#prototypes-for-basic-model-components) describes the prototypes for the basic model components.
Users wishing to write libraries to hold basic model components must furnish functions of these prototypes that perform the basic model component computations.
