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

default: vignettes fresh

vignettes: faq.html mif2.html mif2.R He2010.html He2010.R oaxaca.html oaxaca.R R_v_C.html R_v_C.R upgrade_guide.html C_API.html \
pompjss.pdf getting_started.html getting_started.R

%.Rout: %.R
	R < $*.R > $*.Rout

%.html: %.Rmd
	PATH=/usr/lib/rstudio/bin/pandoc:$$PATH \
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
	$(RM) -r cache figure *.Rout
