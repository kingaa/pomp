RSCRIPT = $(R_HOME)/bin/Rscript --vanilla
REXE = $(R_HOME)/bin/R --vanilla
LATEX = latex
BIBTEX = bibtex
PDFLATEX = pdflatex

default: vignettes clean

vignettes: advanced_topics_in_pomp.pdf advanced_topics_in_pomp.R intro_to_pomp.pdf intro_to_pomp.R

%.R: %.Rnw
	$(RSCRIPT) -e 'library(knitr); purl("$*.Rnw")'
#	$(REXE) CMD Stangle $*.Rnw

%.tex: %.Rnw
	$(RSCRIPT) -e 'library(knitr); knit("$*.Rnw")'
#	$(REXE) CMD Sweave $*.Rnw

%.pdf: %.tex
	$(PDFLATEX) $*
	-$(BIBTEX) $*
	$(PDFLATEX) $*
	$(PDFLATEX) $*
#	$(RSCRIPT) -e "tools::compactPDF(\"$*.pdf\")";

clean:
	$(RM) *.tex *.log *.aux *.blg *.bbl *.out *.Rout *.toc *.lof *.lot *.brf
	$(RM) *.c *.o *.so
	$(RM) Rplots.pdf figure/*
