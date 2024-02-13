---
title: Changes to the way certain functions manipulate the RNG seed
layout: pomp
---

- The `seed` argument to `pfilter` is now deprecated and will be removed soon.
Using it generates a warning.
Equivalent functionality is provided via the new functions `freeze`, `bake`, or `stew`.

- The `seed` argument of `bsmc` and `bsmc2` has been removed.
Its use now generates a warning, stating that it has no effect.
