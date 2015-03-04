RSCRIPT = $(R_HOME)/bin/Rscript --vanilla
REXE = $(R_HOME)/bin/R --vanilla
LATEX = latex
BIBTEX = bibtex
PDFLATEX = pdflatex
RM = rm -f

default: vignettes clean

vignettes: getting_started.html getting_started.R

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
