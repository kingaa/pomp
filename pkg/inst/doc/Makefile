LATEX = latex
PDFLATEX = pdflatex

DVIPS = dvips -Ppdf
EPSTOPDF = epstopdf

default: vignettes clean

vignettes: compiled_code_in_pomp.pdf intro_to_pomp.pdf

%.pdf: %.tex
	$(PDFLATEX) $*
	$(PDFLATEX) $*

%.pdf: %.eps
	$(EPSTOPDF) $*.eps --outfile=$*.pdf

clean:
	$(RM) *.tex *.log *.aux *.out *.Rout *.toc *.lof *.lot *-???.pdf Rplots.ps
