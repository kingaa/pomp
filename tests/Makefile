REXE = R --vanilla
RCMD = $(REXE) CMD

SAVED = $(shell ls *.Rout.save)

default: $(SAVED) clean

%.Rout.save: %.R
	$(RCMD) BATCH --no-timing --no-save $*.R $*.Rout.save
clean:
	$(RM) *.pdf
