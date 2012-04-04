RSCRIPT = $(R_HOME)/bin/Rscript --vanilla
LATEX = latex
BIBTEX = bibtex
PDFLATEX = pdflatex

DVIPS = dvips -Ppdf
EPSTOPDF = epstopdf

default: vignettes compress clean

compress: R_QPDF=qpdf
compress: R_GSCMD=gs
compress: GS_QUALITY=ebook

compress: advanced_topics_in_pomp.pdf intro_to_pomp.pdf manual.pdf
	for i in $^; do \
	 echo "compressing" $$i; \
	 $(RSCRIPT) -e "tools::compactPDF(\"$$i\")"; \
	done 

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
