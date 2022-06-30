# Preparing a CRAN submission

## Checks

- check [CRAN test results](https://cran.r-project.org/web/checks/check_results_pomp.html)
- check [package issues page](https://github.com/kingaa/pomp/issues)
- examine `R-CMD-CHECK` logs on Github Actions
- copy old CRAN branch to new one: `git branch v4.3 v4.2`
- rebase: `git rebase master v4.3`
  - consolidate `NEWS.md`
  - `make NEWS`
- verify that `.Rbuildignore` contains
```
tests
README\.md$
TODO\.md$
^\.gitignore$
^\.travis\.yml$
^codecov\.yml$
^.*\.Rproj$
^\.Rproj\.user$
```
- perform extra checks:
  - `make xxcheck`
  - examine `pomp-Ex.Rout` for `valgrind` output
  - check with R-devel:
```
module switch R/devel
make xcheck
```

## Submitting

- Use the [webform](https://xmpalantir.wu.ac.at/cransubmit/) for submission.
- **Submission note text:**
  The package has been checked under R-4.0.0, R-release, and R-devel on linux, osx, and windows.  I have checked for memory leaks using 'valgrind' and have verified that this release does not break any of the reverse dependencies listed on CRAN. This release fixes the one NOTE mentioned on the CRAN checks page.

## Announcement

- Prepare an announcement for the blog: `www/vignettes/_post/`
- Update manual pages:
  - `make htmlhelp`
  - Then push changes to `kingaa.github.io`
- Make `X.X.0.0` release for Github site.
