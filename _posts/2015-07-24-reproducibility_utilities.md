---
title: Reproducibility utilities
layout: pomp
---

On cooking shows, recipes requiring lengthy baking or stewing are prepared beforehand.
The `bake` and `stew` functions perform analogously, performing an **R** computation and storing the result in a named file.
<!--more-->
If the function is called again and the file is present, the computation is not executed; rather, the results are loaded from the file in which they were previously stored.
Moreover, via their optional `seed` and `kind` arguments, `bake` and `stew` can control the pseudorandom-number generator (RNG) for greater reproducibility.
After the computation is finished, these functions restore the pre-existing RNG state to avoid side effects.

The `freeze` function doesn't save results, but does set the RNG state to the specified value and restore it after the computation is complete.
