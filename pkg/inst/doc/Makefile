LATEX = latex
BIBTEX = bibtex
PDFLATEX = pdflatex

DVIPS = dvips -Ppdf
EPSTOPDF = epstopdf

default: vignettes clean

vignettes: advanced_topics_in_pomp.pdf intro_to_pomp.pdf

%.pdf: %.tex
	$(PDFLATEX) $*
	-$(BIBTEX) $*
	$(PDFLATEX) $*
	$(PDFLATEX) $*

%.pdf: %.eps
	$(EPSTOPDF) $*.eps --outfile=$*.pdf

clean:
	$(RM) *.tex *.log *.aux *.blg *.bbl *.out *.Rout *.toc *.lof *.lot 
	$(RM) Rplots.pdf
	$(RM) intro_to_pomp-*.pdf advanced_topics_in_pomp-*.pdf 
	$(RM) intro_to_pomp-*.png advanced_topics_in_pomp-*.png
