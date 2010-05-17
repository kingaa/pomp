INSTALLDIR = $(HOME)

LATEX = latex
PDFLATEX = pdflatex
BIBTEX = bibtex
MAKEIDX = makeindex
DVIPS = dvips -t letter -Ppdf
RM = rm -f
CP = cp
TAR = tar
TARFLAGS = zcvf
PSTOPDF = ps2pdf
EPSTOPDF = epstopdf

CC = gcc
CXX = g++
CFLAGS = -O2 -I. -Wall -static
LDLIBS = -L. -L/usr/X11/lib 
FC = g77
FFLAGS = -O2 -fno-underscoring
AR = ar
ARFLAGS = rs
INSTALL = install

AWK = gawk
L2H = latex2html
L2HFLAGS = -scalable_fonts -style "whitestyle.css" -short_extn -local_icons -noinfo -split 0 -no_navigation -nosubdir -noaddress -nofootnode
TIDY = tidy
F2D = fig2dev
GZIP = gzip
ZIP = zip

PSFIX = /usr/local/bin/psfix

GRAPH = graph
GRAPHFLAGS = -Tfig -K2 --max-line 256 --title-font-name "times-bolditalic" -F"times-bold" -W0.005 -B

TABLES = 

FIGS = 

EPSFIGS = $(addsuffix .eps, $(FIGS))

PDFFIGS = $(addsuffix .pdf, $(FIGS))

%.ps: %.dvi
	$(DVIPS) $* -o $*.ps

%.dvi: %.tex
	$(LATEX) $*
	-$(BIBTEX) $*
	$(LATEX) $*
	$(LATEX) $*

%.pdf: %.tex
	$(PDFLATEX) $*
	-$(BIBTEX) $*
	$(PDFLATEX) $*
	$(PDFLATEX) $*

%.pdf: %.ps
	$(PSTOPDF) $*.ps $*.pdf

%.idx: %.tex
	-$(LATEX) $*

%.ind: %.idx
	$(MAKEIDX) $*

%.pdf: %.eps
	$(EPSTOPDF) $*.eps --outfile=$*.pdf

%.pdf: %.fig
	$(F2D) -Lpdf $*.fig $@

%.eps: %.fig
	$(F2D) -Leps $*.fig $@

%.gif: %.fig
	$(F2D) -Lgif -t"#ffffff" $*.fig $@

%.png: %.fig
	$(F2D) -Lpng $*.fig $@

%.pdf: %.fig
	$(F2D) -Lpdf $*.fig $@

%.pcx: %.fig
	$(F2D) -Lpcx $*.fig $@

%.htm: %.tex %.pdf
	$(L2H) $(L2HFLAGS) $*
	$(AWK) '{gsub("H2", "H3"); print;}' $*.htm > $*_1.htm
	$(AWK) '{gsub("<BR><HR>", ""); print;}' $*_1.htm > $*_2.htm
	-$(TIDY) -i -w --show-body-only yes $*_2.htm > $*.htm
	$(RM) $*_1.htm $*_2.htm

%.ps: %.nb.ps
	$(PSFIX) $*.nb.ps > $@

%.gz: %
	$(GZIP) < $* > $*.gz

%.tgz: %
	$(TAR) $(TARFLAGS) $*.tgz $*

%: %.c
	$(CC) $(CFLAGS) $*.c $(LDLIBS) -o $*

%: %.cc
	$(CXX) $(CFLAGS) $*.cc $(LDLIBS) -o $*

%: %.f
	$(FC) $(FFLAGS) $*.f $(LDLIBS) -o $*

(%.o): %.c
	$(CC) -c $(CFLAGS) $*.c -o $*.o
	$(AR) $(ARFLAGS) $@ $*.o
	$(RM) $*.o 

(%.o): %.cc
	$(CXX) -c $(CFLAGS) $*.cc -o $*.o
	$(AR) $(ARFLAGS) $@ $*.o
	$(RM) $*.o 

(%.o): %.f
	$(FC) -c $(FFLAGS) $*.f -o $*.o
	$(AR) $(ARFLAGS) $@ $*.o
	$(RM) $*.o

$(LIB):	$(LIB)($(addsuffix .o, $(basename $(SRC))))

clean:
	$(RM) *.log *.blg *.ilg *.aux *.lof *.lot *.toc *.idx *.ttt *.fff *.out internals.pl labels.pl

fresh: clean
	$(RM) *.bbl *.ind *.dvi *.ps
