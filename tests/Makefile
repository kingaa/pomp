REXE = R

SAVED = $(wildcard *.Rout.save *.Rout)

default: $(SAVED) clean

%.Rout.save: %.R
	$(REXE) < $*.R > $*.Rout.save 2>&1

%.Rout: %.R
	$(REXE) < $*.R > $*.Rout 2>&1

clean:
	$(RM) *.c *.o *.so 

fresh: clean
	$(RM) *.png
