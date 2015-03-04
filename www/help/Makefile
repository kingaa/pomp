RSCRIPT = $(R_HOME)/bin/Rscript --vanilla
REXE = $(R_HOME)/bin/R --vanilla --slave
RM = rm -f

default: clean help

help:
	$(REXE) -e 'knitr::knit_rd("pomp")'

clean:
	$(RM) *.html
