----------------

# Instructions for building, checking, and installing

## Makefile usage

A `Makefile` is provided.  Most of the rules are defined in the `rules.mk` file.
The following key variables are defined in the `Makefile`.

- `REPODIR`: the name of the directory in which a local copy of the Git repo used for publishing source and documentation is located.
- `INCLUDES`: names of header files that will be copied to `inst/includes`
- `HEADERS`: names of header files that must be rebuilt whenever source code changes.

The key targets are:

- `roxy`: runs `devtools::document()` to build the help pages (in `man/`), the package `NAMESPACE`, and the collation order in `DESCRIPTION`.
- `dist`: builds the package source tarball.
- `install`: installs the package locally for testing.
- `tests`: runs `make` in the `tests/` directory, causing tests to be updated as needed.
- `clean`: cleans the directory of files created during builds and tests.
- `fresh`: resets the directory for a fresh rebuild.

- `htmlhelp`: builds the package manual in HTML and PDF formats.
  This command also causes `doxygen` to be run, which produces HTML documentation of the C/C++ source code.
  These documentation files are installed into the `REPODIR/manuals` directory.
- `www`: installs the package and runs `make` in the `www` directory.
  This causes vignettes to be updated, among other things.
- `NEWS`: builds the plain-text package NEWS file.
- `instdocs`: runs `make` in the `inst/docs` directory.

- `check`: runs `devtools::check()`.
- `qcheck`: like `check`, but tests in `tests/` are not run.
- `qqcheck`: like `qcheck`, but checks for code-documentation mismatch and examples are not run.
- `xcheck`: runs `devtools::check(cran=TRUE)`, i.e., the additional checks for CRAN suitability are also run.
- `ycheck`: runs `xcheck` but tests and examples that are ordinarily not run (because expensive) are run.
- `vcheck`: runs `check` and then tests for memory leaks using `valgrind`.
  One must examine the output (printed to the console and also stored in a file, the name of which ends in `-Ex.Rout`.
- `revdeps`: causes checks of reverse-dependent packages to run.
  The outputs are stored under the `revdep` directory.
- `rchk`: runs Kalibera's `rchk` utility.
  This uses `docker` to download and run a Docker container.
  The stdout for this is stored in `rchk.out`.

- `session`: installs the package locally and then runs an **R** session for interactive testing.
  This session is controlled by the `RSESSION` environment variable.
  By default, this runs an `emacs -f R` session.
- `rsession`: like `session`, but a naked **R** session is started: `RSESSION=R`.
- `debug`: like `rsession`, but the **R** session is started under a debugger.
  By default, the `gdb` debugger is used: `RSESSION=R -d gdb`.

- `publish`: causes source and binary tarballs to be installed in the REPODIR respository.
- `covr`: runs `covr::package_coverage()` to evaluate unit-test coverage.
  Results are stored in `covr.rds`.
- `vcovr`: like `covr`, but results are displayed in a browser. 

----------------

## Github Actions

Several Github Actions are defined in `.github/workflows`:

- `r-cmd-check`: runs checks on various platforms.
- `test-coverage`: runs `covr::codecov` to upload unit-test coverage information to [codecov.io](https://codecov.io).
- `binary-build`: causes binary tarballs to be built for OS X and Windows platforms.

These must be updated manually from time to time as the Actions on which they depend change.

----------------
