ifneq ($(R_HOME),)
	RSCRIPT = $(R_HOME)/bin/Rscript --vanilla
	REXE = $(R_HOME)/bin/R --vanilla --slave
else
	RSCRIPT = Rscript --vanilla
	REXE = R --vanilla --slave
endif
LATEX = latex
BIBTEX = bibtex
PDFLATEX = pdflatex
RM = rm -f

default: vignettes clean

vignettes: FAQ.html He2010.R He2010.html oaxaca.R oaxaca.html \
R_v_C.R R_v_C.html upgradeGuide.html C_API.html \
pompjss.R pompjss.pdf getting_started.R getting_started.html

%.Rout: %.R
	R < $*.R > $*.Rout

%.html: %.Rmd
	$(RSCRIPT) -e "rmarkdown::render(\"$*.Rmd\")"

%.R: %.Rnw
	$(RSCRIPT) -e 'library(knitr); purl("$*.Rnw")'

%.tex: %.Rnw
	$(RSCRIPT) -e 'library(knitr); knit("$*.Rnw")'

%.R: %.Rmd
	$(RSCRIPT) -e 'library(knitr); purl("$*.Rmd")'

%.pdf: %.tex
	$(PDFLATEX) $*
	-$(BIBTEX) $*
	$(PDFLATEX) $*
	$(PDFLATEX) $*

clean:
	$(RM) *.tex *.log *.aux *.blg *.bbl *.out *.toc *.lof *.lot *.brf
	$(RM) *.c *.o *.so *.bak
	$(RM) Rplots.pdf

fresh: clean
	$(RM) -r tmp *.Rout
