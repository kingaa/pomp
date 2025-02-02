---
layout: pomp
title: pomp version 6.1 released
---

**pomp** version 6.1 has been released to CRAN and will soon appear on [a mirror near you](https://cran.r-project.org/mirrors.html).
This release contains breaking changes as well as some additional features.

### User-visible changes

#### New interface for `userdata` becomes mandatory

Since version 5.8.4, the manner in which one provides extra elements to basic model components (i.e., beyond parameters, latent state variables, covariates, time, and observations: the so-called [`userdata`](https://kingaa.github.io/manuals/pomp/html/userdata.html)) has changed. 
During a grace period, the old method still worked, though it generated a warning.
In versions 6.X, an error will be generated.
To supply additional elements to the basic model components, pass them as elements of a named list via the `userdata` argument, which can be furnished to any [elementary algorithm](https://kingaa.github.io/manuals/pomp/html/elementary_algorithms.html) or [estimation algorithm](https://kingaa.github.io/manuals/pomp/html/estimation_algorithms.html), and of course, to the [pomp constructor](https://kingaa.github.io/manuals/pomp/html/pomp.html) itself.

#### Passing arguments by position now results in an error in most cases

In calls to **pomp** elementary and inference algorithms, it is now necessary to pass arguments *by name* and not by position.
This has always been good practice, but from this release, calls that rely on the position of arguments will typically generate errors.

### Feature enhancements

#### Initial value parameters in `mif2`

It is now possible to specify more than one lag in the `ivp` function, which is evaluated only when the `mif2` perturbations are specified.
See `?mif2`.

#### Keeping a database of explorations

In conducting an extensive exploration of a likelihood surface, it is useful to maintain a database of places visited, together with associated likelihood values.
The new function [`append_data`](https://kingaa.github.io/manuals/pomp/html/bake.html#:~:text=append_data) assists in this.
It simply appends a given data frame to an existing CSV file (creating the file if it does not exist).

#### Expectation of an Euler-multinomial random variable

The new function [`eeulermultinom`](https://kingaa.github.io/manuals/pomp/html/eulermultinom.html#:~:text=eeulermultinom) gives the expectation of an Euler-multinomial random variable.
There is also an interface to this function in the [C API](C_API.html#expectation-of-an-euler-multinomial-random-variable).

### Other changes

The `save.states` option to `pfilter` has changed.
See [`pfilter`](https://kingaa.github.io/manuals/pomp/html/pfilter.html) and [`saved_states`](https://kingaa.github.io/manuals/pomp/html/saved_states.html) for details.
The deprecated options will still work for the present, but will generate a warning, with advice about how to change.
Ultimately, these options will be removed.
