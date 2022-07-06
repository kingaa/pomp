---
date: 8 July 2022
layout: pomp
title: version 4.3 released to CRAN
---
    
Version 4.3 of **pomp** is now on CRAN and coming soon to [a mirror near you](https://cran.r-project.org/mirrors.html).
This release contains a few changes to the user interface.
- The archiving functions `bake` and `stew` now use a slightly less exacting comparison of the expression, `expr`, they are furnished.
  **NB:** Running `bake` or `stew` with archives created by earlier versions may result in recomputation.
- Also, `bake` and `stew` now take the argument `dir`, which is the directory holding the archive files.
  By default, this is the current working directory or the value of the global option `pomp.archive.dir`.
- All workhorse functions except `partrans` have new default arguments.
- The `time` method has been extended to `pompList` and related objects.
  **NB:** The behavior of this function may be streamlined in the near future.


In addition, there are some small improvements to the [manual](https://kingaa.github.io/manuals/pomp/).
