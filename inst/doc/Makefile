LATEX = latex
PDFLATEX = pdflatex

DVIPS = dvips -Ppdf
EPSTOPDF = epstopdf

VIGNETTES = random_walk_example.pdf

default: vignettes clean

vignettes: $(VIGNETTES)

%.pdf: %.tex
	$(PDFLATEX) $*
	$(PDFLATEX) $*

%.pdf: %.eps
	$(EPSTOPDF) $*.eps --outfile=$*.pdf

clean:
	$(RM) *.tex *.log *.aux *.out *.Rout *-???.pdf Rplots.ps
