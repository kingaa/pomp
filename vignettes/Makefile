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

vignettes: getting_started.html getting_started.R pompjss.pdf pompjss.R mif2.html mif2.R He2010.html He2010.R nlf.html oaxaca.html oaxaca.R R_v_C.html R_v_C.R upgrade_guide.html

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
#	$(RSCRIPT) -e "tools::compactPDF(\"$*.pdf\")";

clean:
	$(RM) *.tex *.log *.aux *.blg *.bbl *.out *.Rout *.toc *.lof *.lot *.brf
	$(RM) -r cache figure
	$(RM) *.c *.o *.so *.bak
	$(RM) Rplots.pdf

fresh: clean
