REXE = R --vanilla
RCMD = $(REXE) CMD

SAVED = $(shell ls *.Rout.save *.Rout)
DEMOS = $(shell ls ../demo/*R)

default: $(SAVED) clean

demos.Rout: $(DEMOS)

%.Rout.save: %.R
	$(RCMD) BATCH --no-timing --no-save $*.R $*.Rout.save

%.Rout: %.R
	$(RCMD) BATCH --no-timing --no-save $*.R $*.Rout

clean:
	$(RM) *.c *.o *.so 

fresh: clean
	$(RM) *.png
