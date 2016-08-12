REXE = R --vanilla
RCMD = $(REXE) CMD

SAVED = $(shell ls *.Rout.save)

default: $(SAVED) demos.Rout verbosity.Rout link.Rout clean

%.Rout.save: %.R
	$(RCMD) BATCH --no-timing --no-save $*.R $*.Rout.save

%.Rout: %.R
	$(RCMD) BATCH --no-timing --no-save $*.R $*.Rout

clean:
	$(RM) *.pdf *.Rout *.c *.o *.so 

fresh: clean
	$(RM) *.png
