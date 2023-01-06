REXE = R -s --vanilla
REXE_ALT = R -s
RSESSION = emacs -f R
RCMD = $(REXE) CMD
RCMD_ALT = R -s --no-save --no-restore CMD
RSCRIPT = Rscript --vanilla
REPODIR = ../www
MANUALDIR = ../www/manuals/$(PKG)

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
SOURCE = $(sort $(wildcard R/*R src/*.c src/*.h data/* examples/*))
CSOURCE = $(sort $(wildcard src/*.c))
TESTS = $(sort $(wildcard tests/*R))
INSTDOCS = $(sort $(wildcard inst/doc/*))
SESSION_PKGS = datasets,utils,grDevices,graphics,stats,methods,tidyverse,$(PKG)

default:
	@echo $(PKGVERS)

.PHONY: binary check clean covr debug default fresh \
htmlhelp manual news publish qcheck qqcheck \
revdeps rhub rsession session vignettes win wind xcheck \
xcovr xxcheck ycheck

.dist manual vignettes: export R_QPDF=qpdf
.headers: export LC_COLLATE=C
.roxy .headers .dist manual vignettes: export R_HOME=$(shell $(REXE) RHOME)
check xcheck xxcheck: export FULL_TESTS=yes
.dist .tests session check xcheck xxcheck: export R_KEEP_PKG_SOURCE=yes
revdeps .tests xcheck: export R_PROFILE_USER=$(CURDIR)/.Rprofile
.tests session xxcheck vignettes data manual: export R_LIBS=$(CURDIR)/library
xcheck: export _R_CHECK_DEPENDS_ONLY_=true
debug: export RSESSION=R -d gdb
rsession: export RSESSION=R

inst/include/%.h: src/%.h
	$(CP) $^ $@

dist: .dist

install: .install

instdocs: .instdocs

roxy: .roxy

htmlhelp: install news manual
	rsync --delete -a library/$(PKG)/html/ $(MANUALDIR)/html
	rsync --delete --exclude=aliases.rds --exclude=paths.rds --exclude=$(PKG).rdb --exclude=$(PKG).rdx --exclude=macros -a library/$(PKG)/help/ $(MANUALDIR)/help
	(cd $(MANUALDIR); (cat links.ed && echo w ) | ed - html/00Index.html)
	$(CP) $(PKG).pdf $(MANUALDIR)
	$(CP) $(REPODIR)/assets/R.css $(MANUALDIR)/html

vignettes: manual install
	$(MAKE)	-C www

news: library/$(PKG)/html/NEWS.html

NEWS: .NEWS

.instdocs: $(INSTDOCS)
	$(MAKE) -C inst/doc
	$(TOUCH) $@

.NEWS: inst/NEWS
	$(TOUCH) $@

inst/NEWS: inst/NEWS.Rd
	$(RCMD) Rdconv -t txt $^ -o $@

library/$(PKG)/html/NEWS.html: inst/NEWS.Rd
	$(RCMD) Rdconv -t html $^ -o $@

session: install
	export R_DEFAULT_PACKAGES=$(SESSION_PKGS) && \
	exec $(RSESSION)

debug: session

rsession: session

revdeps: .dist
	mkdir -p revdep
	$(CP) $(PKGVERS).tar.gz revdep
	$(REXE_ALT) -e "tools::check_packages_in_dir(\"revdep\",check_args=\"--as-cran\",reverse=list(which=\"most\"))"

.roxy: .source .headers
	$(REXE) -e "devtools::document()"
	$(TOUCH) $@

.headers: $(HEADERS)
	$(TOUCH) $@

.includes: $(INCLUDES)
	$(TOUCH) $@

.source: $(SOURCE)
	$(TOUCH) $@

.testsource: $(TESTS)
	$(TOUCH) $@

.dist: .roxy .NEWS .instdocs .source .testsource .includes .headers
	$(RCMD) build --force --no-manual --resave-data --compact-vignettes=both --md5 .
	$(TOUCH) $@

binary: dist
	mkdir -p plib
	$(RCMD) INSTALL --build --library=plib --preclean --clean $(PKGVERS).tar.gz
	rm -rf plib

publish: dist manual htmlhelp
	$(RSCRIPT) -e 'drat::insertPackage("$(PKGVERS).tar.gz",repodir="$(REPODIR)",action="prune")'
	-$(RSCRIPT) -e 'drat::insertPackage("$(PKGVERS).tgz",repodir="$(REPODIR)",action="prune")'
	-$(RSCRIPT) -e 'drat::insertPackage("$(PKGVERS).zip",repodir="$(REPODIR)",action="prune")'

rhub:
	$(REXE) -e 'library(rhub); check_for_cran(); check_on_windows(); check(platform="macos-highsierra-release-cran");'

covr: covr.rds

covr.rds: DESCRIPTION
	$(REXE) -e 'library(covr); package_coverage(type="all") -> cov; report(cov,file="covr.html",browse=TRUE); saveRDS(cov,file="covr.rds")'

xcovr: covr
	$(REXE) -e 'library(covr); readRDS("covr.rds") -> cov; codecov(coverage=cov,quiet=FALSE)'

win: dist
	curl -T $(PKGVERS).tar.gz ftp://win-builder.r-project.org/R-release/

wind: dist
	curl -T $(PKGVERS).tar.gz ftp://win-builder.r-project.org/R-devel/

check: dist
	mkdir -p check
	$(RCMD) check --no-stop-on-test-error --library=check -o check $(PKGVERS).tar.gz

qcheck: dist
	mkdir -p check
	$(RCMD) check --library=check -o check --no-vignettes --no-tests $(PKGVERS).tar.gz

qqcheck: dist
	mkdir -p check
	$(RCMD) check --library=check -o check --no-codoc --no-examples --no-vignettes --no-manual --no-tests $(PKGVERS).tar.gz

xcheck: dist
	mkdir -p check library
	$(RCMD_ALT) check --no-stop-on-test-error --as-cran --library=library -o check $(PKGVERS).tar.gz

xxcheck: install xcheck
	mkdir -p check
	$(REXE) -d "valgrind --tool=memcheck --track-origins=yes --leak-check=full" < check/$(PKG).Rcheck/$(PKG)-Ex.R 2>&1 | tee $(PKG)-Ex.Rout

ycheck: dist install
	mkdir -p check
	$(RCMD_ALT) check --run-dontrun --run-donttest --as-cran --library=library -o check $(PKGVERS).tar.gz

manual: install $(PKG).pdf

$(PKG).pdf: $(SOURCE)
	$(RCMD) Rd2pdf --internals --no-description --no-preview --pdf --force -o $(PKG).pdf .
	$(RSCRIPT) -e "tools::compactPDF(\"$(PKG).pdf\")";

tests: .tests

.tests: .install .testsource
	$(MAKE) -C tests
	$(TOUCH) $@

install: .install

.install: .roxy .NEWS .instdocs .source .includes .headers
	mkdir -p library
	$(RCMD) INSTALL --html --library=library .
	$(TOUCH) .install

%.tex: %.Rnw
	$(RSCRIPT) -e "library(knitr); knit(\"$*.Rnw\")"

%.R: %.Rnw
	$(RSCRIPT) -e "library(knitr); purl(\"$*.Rnw\")"

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
	PATH=/usr/lib/rstudio/bin/pandoc:$$PATH \
	Rscript --vanilla -e "rmarkdown::render(\"$*.Rmd\")"

%.html: %.md
	PATH=/usr/lib/rstudio/bin/pandoc:$$PATH \
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
	$(RM) .dist

fresh: clean
	$(RM) .headers .includes .NEWS .instdocs
	$(RM) .install .roxy .source .testsource .roxy .tests
	$(RM) -r library
	$(RM) -r revdep/*.Rcheck revdep/Library revdep/*.tar.gz
