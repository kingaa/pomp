# Preparing a CRAN submission

## Checks

- check [CRAN test results](https://cran.r-project.org/web/checks/check_results_pomp.html)
- check [package issues page](https://github.com/kingaa/pomp/issues)
- examine `R-CMD-CHECK` logs on Github Actions
- examine tests coverage on [codecov.com](https://codecov.io/gh/kingaa/pomp).
- copy old CRAN branch to new one and update timestamp:

```
git branch v4.6 v4.5
git rebase master v4.6
git commit --amend --date=now --no-edit
git push github v4.6
```

- rebase: `git rebase master v4.6`
  - consolidate `NEWS.md`
  - `make NEWS`
- verify that `.Rbuildignore` contains
  - `tests`
  - `README\.md$`
  - `TODO\.md$`
  - `CRAN\.md$`
  - `^\.gitignore$`
  - `^\.travis\.yml$`
  - `^codecov\.yml$`
- perform extra checks:
  - `make xxcheck`
  - examine `pomp-Ex.Rout` for `valgrind` output

## Submitting

- Use the [webform](https://xmpalantir.wu.ac.at/cransubmit/) for submission.
- **Submission note text:**  
  The package has been checked under R-4.1.0, R-release, and R-devel on linux, osx, and windows.
  I have checked for memory leaks using 'valgrind'.
  There are two new NOTEs generated, concerning S3 generic/method consistency.
  These are false positives, since the functions in question `filter.mean` and `filter.traj` are themselves S4 generics.
  I have checked that the S3 method dispatch system works properly.
  In any case, these functions have been deprecated as of this version in favor of functions named `filter_mean` and `filter_traj`, so this false positive will go away when these functions are removed in a future release.
  I have verified that this release does not break the reverse dependencies listed on CRAN, with the exception of **spatPomp**.
  I am one of the developers of that package and we have prepared a CRAN release of this package that corrects this problem.
  We will submit this as soon as this version of **pomp** is published to CRAN.

## Announcement

- Prepare an announcement for the blog: `www/vignettes/_post/`
- Update manual pages:
  - `make publish`
  - Then push changes to `kingaa.github.io`
- Make `X.X.0.0` release for Github site.
