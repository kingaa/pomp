RSCRIPT = $(R_HOME)/bin/Rscript --vanilla
REXE = $(R_HOME)/bin/R --vanilla
LATEX = latex
BIBTEX = bibtex
PDFLATEX = pdflatex

default: vignettes clean

vignettes: advanced_topics_in_pomp.pdf advanced_topics_in_pomp.R intro_to_pomp.pdf intro_to_pomp.R

%.R: %.Rnw
	$(REXE) CMD Stangle $*.Rnw

%.tex: %.Rnw
	$(REXE) CMD Sweave $*.Rnw

%.pdf: %.tex
	$(PDFLATEX) $*
	-$(BIBTEX) $*
	$(PDFLATEX) $*
	$(PDFLATEX) $*
	$(RSCRIPT) -e "tools::compactPDF(\"$*.pdf\")";

clean:
	$(RM) *.tex *.log *.aux *.blg *.bbl *.out *.Rout *.toc *.lof *.lot 
	$(RM) *.c *.o *.so
	$(RM) Rplots.pdf
	$(RM) intro_to_pomp-*.pdf advanced_topics_in_pomp-*.pdf 
	$(RM) intro_to_pomp-*.png advanced_topics_in_pomp-*.png
