---
date: 03 March 2023
layout: pomp
title: version 4.7 released to CRAN
---
    
Version 4.7 of **pomp** is now on CRAN;
download it from [a mirror near you](https://cran.r-project.org/mirrors.html).
In this release, some deprecated functions have become defunct.

### Changes in function names

As of version 4.6, all **pomp** functions with names.that.contain.dots were deprecated in favor of functions in snake_case.
That is, every function that had a dot (`.`) in its name was replaced by a function where every dot is replaced by an underscore (`_`).
These functions have now become defunct: using them will result in an error message that tells you what the replacement function is.

To help you adapt your code to the new naming convention, you can download and run the [to_snake_case.R](https://kingaa.github.io/scripts/to_snake_case.R) script.
Its usage is straightforward:  
1. Make a directory and copy all files that you wish to edit into it.
2. In an **R** session, source this script.
3. Call the `to_snake_case()` function with the path to your new directory as its sole argument
4. Examine the differences between the files for correctness.
5. Move the new files back into place.

### Other user-visible changes

- The `plot` function now works when the data contain missing (`NA`) values.
