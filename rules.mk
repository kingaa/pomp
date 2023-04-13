REXE = $(shell which R) -s
RCMD = $(REXE) CMD
MANUALDIR = $(REPODIR)/manuals/$(PKG)

PDFLATEX = pdflatex
BIBTEX = bibtex
MAKEIDX = makeindex

RM = rm -f
CP = cp
TOUCH = touch
INSTALL = install

PKG = $(shell perl -ne 'print $$1 if /Package:\s+((\w+[-\.]?)+)/;' DESCRIPTION)
VERSION = $(shell perl -ne 'print $$1 if /Version:\s+((\d+[-\.]?)+)/;' DESCRIPTION)
PKGVERS = $(PKG)_$(VERSION)
TARBALL = $(PKGVERS).tar.gz
SOURCE = $(sort $(wildcard R/*R src/*.c src/*.h data/* examples/*))
CSOURCE = $(sort $(wildcard src/*.c))
TESTS = $(sort $(wildcard tests/*R))
INSTDOCS = $(sort $(wildcard inst/doc/*))
SESSION_PKGS = datasets,utils,grDevices,graphics,stats,methods,tidyverse,$(PKG)

.PHONY: .check check clean covr debug default fresh \
htmlhelp manual publish qcheck qqcheck \
revdeps rhub rsession session www win wind xcheck \
xcovr vcheck ycheck

.dist manual www: export R_QPDF=qpdf
.headers: export LC_COLLATE=C
.roxy .headers .dist manual www: export R_HOME=$(shell $(REXE) RHOME)
.check: export FULL_TESTS=yes
.dist .tests .session .check: export R_KEEP_PKG_SOURCE=yes
revdeps .session .tests .check: export R_PROFILE_USER=$(CURDIR)/.Rprofile
.tests .session vcheck www manual: export R_LIBS=$(CURDIR)/library
.check: export R_CHECK_ENVIRON=$(CURDIR)/tools/check.env
session: RSESSION = emacs -f R
debug: RSESSION = R -d gdb
rsession: RSESSION = R

default: .roxy .NEWS .instdocs .source .includes .headers
	@echo $(PKGVERS)

roxy: .roxy

dist: .dist

install: .install

instdocs: .instdocs

COMMON_CHECK_ARGS = document=FALSE,vignettes=FALSE,clean_doc=FALSE,check_dir="check"

check: CHECK = devtools::check($(COMMON_CHECK_ARGS),cran=FALSE)
qcheck: CHECK = devtools::check($(COMMON_CHECK_ARGS),cran=FALSE,\
args=c("--no-tests"))
qqcheck: CHECK = devtools::check($(COMMON_CHECK_ARGS),cran=FALSE,\
args=c("--no-tests","--no-codoc","--no-examples"))
xcheck: CHECK = devtools::check($(COMMON_CHECK_ARGS),cran=TRUE,\
env_vars=c("_R_CHECK_DEPENDS_ONLY_"="TRUE"))
ycheck: CHECK = devtools::check($(COMMON_CHECK_ARGS),cran=TRUE,\
args=c("--run-dontrun","--run-donttest"))

INSTALLCMD = devtools::install(args=c("--preclean","--html","--library=library"))

check xcheck ycheck qcheck qqcheck: .check

vcheck: check/$(PKG).Rcheck/$(PKG)-Ex.R
	$(REXE) -d "valgrind -s --tool=memcheck --track-origins=yes --leak-check=full"\
	< $^ 2>&1 | tee $(PKG)-Ex.Rout

NEWS: .NEWS

.instdocs: $(INSTDOCS)
	$(MAKE) -C inst/doc
	$(TOUCH) $@

.NEWS: inst/NEWS
	$(TOUCH) $@

.headers: $(HEADERS)
	$(TOUCH) $@

.includes: $(INCLUDES)
	$(TOUCH) $@

.source: $(SOURCE)
	$(TOUCH) $@

.testsource: $(TESTS)
	$(TOUCH) $@

.roxy: .source .headers
	$(REXE) -e "devtools::document()"
	$(TOUCH) $@

.check: .roxy .NEWS .instdocs .includes
	$(REXE) -e '$(CHECK)'

.dist: .roxy .NEWS .instdocs .testsource .includes
	$(RCMD) build --force --no-manual --resave-data --compact-vignettes=both --md5 .
	$(TOUCH) $@

.install: .roxy .NEWS .instdocs .source .includes .headers
	mkdir -p library
	$(REXE) -e '$(INSTALLCMD)'
	$(RCMD) Rdconv -t html inst/NEWS.Rd -o library/$(PKG)/html/NEWS.html
	$(TOUCH) .install

.session: .install
	export R_DEFAULT_PACKAGES=$(SESSION_PKGS) && \
	exec $(RSESSION)

inst/NEWS: inst/NEWS.Rd
	$(RCMD) Rdconv -t txt $^ -o $@

htmlhelp: install manual
	rsync --delete -a library/$(PKG)/html/ $(MANUALDIR)/html
	rsync --delete --exclude=aliases.rds --exclude=paths.rds --exclude=$(PKG).rdb --exclude=$(PKG).rdx --exclude=macros -a library/$(PKG)/help/ $(MANUALDIR)/help
	(cd $(MANUALDIR); (cat links.ed && echo w ) | ed - html/00Index.html)
	$(CP) $(PKG).pdf $(MANUALDIR)
	$(CP) $(REPODIR)/assets/R.css $(MANUALDIR)/html

www: install
	$(MAKE)	-C www

session debug rsession: .session

revdeps: .dist
	mkdir -p revdep
	$(CP) $(TARBALL) revdep
	$(REXE) -e "tools::check_packages_in_dir(\"revdep\",check_args=\"--as-cran\",reverse=list(which=\"most\"))"

publish: dist manual htmlhelp
	$(REXE) -e 'drat::insertPackage("$(PKGVERS).tar.gz",repodir="$(REPODIR)",action="prune")'
	-$(REXE) -e 'drat::insertPackage("$(PKGVERS).tgz",repodir="$(REPODIR)",action="prune")'
	-$(REXE) -e 'drat::insertPackage("$(PKGVERS).zip",repodir="$(REPODIR)",action="prune")'

rhub:
	$(REXE) -e 'library(rhub); check_for_cran(); check_on_windows(); check(platform="macos-highsierra-release-cran");'

covr: covr.rds

covr.rds: DESCRIPTION
	$(REXE) -e 'library(covr); package_coverage(type="all") -> cov; report(cov,file="covr.html",browse=TRUE); saveRDS(cov,file="covr.rds")'

xcovr: covr
	$(REXE) -e 'library(covr); readRDS("covr.rds") -> cov; codecov(coverage=cov,quiet=FALSE)'

win: dist
	curl -T $(TARBALL) ftp://win-builder.r-project.org/R-release/

wind: dist
	curl -T $(TARBALL) ftp://win-builder.r-project.org/R-devel/

manual: install $(PKG).pdf

$(PKG).pdf: $(SOURCE)
	$(RCMD) Rd2pdf --internals --no-description --no-preview --pdf --force -o $(PKG).pdf .
	$(REXE) -e "tools::compactPDF(\"$(PKG).pdf\")";

tests: .tests

.tests: .install .testsource
	$(MAKE) -j10 -C tests
	$(TOUCH) $@

install: .install

inst/include/%.h: src/%.h
	$(CP) $^ $@

%.tex: %.Rnw
	$(REXE) -e "library(knitr); knit(\"$*.Rnw\")"

%.R: %.Rnw
	$(REXE) -e "library(knitr); purl(\"$*.Rnw\")"

%.pdf: %.tex
	$(PDFLATEX) $*
	-$(BIBTEX) $*
	$(PDFLATEX) $*
	$(PDFLATEX) $*

%.bbl: %.tex
	-$(PDFLATEX) $*
	$(BIBTEX) $*

%.idx: %.tex
	-$(PDFLATEX) $*

%.ind: %.idx
	$(MAKEIDX) $*

%.html: %.Rmd
	Rscript --vanilla -e "rmarkdown::render(\"$*.Rmd\")"

%.html: %.md
	Rscript --vanilla -e "rmarkdown::render(\"$*.md\")"

%.R: %.Rmd
	Rscript --vanilla -e "knitr::purl(\"$*.Rmd\",output=\"$*.R\",documentation=2)"

clean:
	$(RM) -r check
	$(RM) src/*.o src/*.so src/symbols.rds
	$(RM) -r lib
	$(RM) -r *-Ex.Rout *-Ex.timings *-Ex.pdf Rplots.*
	$(RM) *.tar.gz $(PKGVERS).zip $(PKGVERS).tgz $(PKG).pdf
	$(MAKE) -C www clean
	$(MAKE) -C inst/doc clean
	$(MAKE) -C tests clean
	$(MAKE) -C revdep clean
	$(RM) .dist

fresh: clean
	$(RM) .headers .includes .NEWS .instdocs
	$(RM) .install .roxy .source .testsource .roxy .tests
	$(RM) -r library
