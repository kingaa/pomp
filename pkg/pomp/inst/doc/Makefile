RSCRIPT = $(R_HOME)/bin/Rscript --vanilla
REXE = $(R_HOME)/bin/R --vanilla
LATEX = latex
BIBTEX = bibtex
PDFLATEX = pdflatex

default: vignettes clean

vignettes: advanced_topics_in_pomp.pdf intro_to_pomp.pdf

%.tex: %.Rnw
	$(REXE) CMD Sweave $*

%.pdf: %.tex
	$(PDFLATEX) $*
	-$(BIBTEX) $*
	$(PDFLATEX) $*
	$(PDFLATEX) $*
	$(RSCRIPT) -e "tools::compactPDF(\"$*.pdf\")";

clean:
	$(RM) *.tex *.log *.aux *.blg *.bbl *.out *.Rout *.toc *.lof *.lot 
	$(RM) Rplots.pdf
	$(RM) intro_to_pomp-*.pdf advanced_topics_in_pomp-*.pdf 
	$(RM) intro_to_pomp-*.png advanced_topics_in_pomp-*.png
