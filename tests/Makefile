REXE = R
RCMD = $(REXE) CMD

SAVED = $(shell ls *.Rout.save *.Rout)
DEMOS = $(shell ls ../demo/*R)

default: $(SAVED) clean

demos.Rout: $(DEMOS)

%.Rout.save: %.R
	$(REXE) < $*.R > $*.Rout.save 2>&1

%.Rout: %.R
	$(REXE) < $*.R > $*.Rout 2>&1

clean:
	$(RM) *.c *.o *.so 

fresh: clean
	$(RM) *.png
