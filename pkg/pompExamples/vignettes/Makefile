RSCRIPT = $(R_HOME)/bin/Rscript --vanilla
REXE = $(R_HOME)/bin/R --vanilla
LATEX = latex
BIBTEX = bibtex
PDFLATEX = pdflatex

default: vignettes clean

vignettes: budmoth-model.pdf pertussis-model.pdf

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
	$(RM) budmoth-model-*.pdf pertussis-model-*.pdf 
	$(RM) budmoth-model-*.png pertussis-model-*.png
