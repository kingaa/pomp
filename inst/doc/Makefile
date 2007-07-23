LATEX = latex
PDFLATEX = pdflatex

DVIPS = dvips -Ppdf
EPSTOPDF = epstopdf

default: vignettes clean
	mv ou2.rda ../../data

vignettes: random_walk.pdf

%.pdf: %.tex
	$(PDFLATEX) $*
	$(PDFLATEX) $*

%.pdf: %.eps
	$(EPSTOPDF) $*.eps --outfile=$*.pdf

clean:
	$(RM) *.tex *.log *.aux *.out *.Rout *-???.pdf Rplots.ps
