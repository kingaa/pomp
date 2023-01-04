REVDEPS=spaero epimdr epimdr2 CollocInfer spatPomp DTAT
INCLUDES=inst/include/pomp.h inst/include/pomp_defines.h
HEADERS=src/pomp_decls.h

include rules.mk

src/pomp_decls.h: $(CSOURCE)
	cproto -I $(R_HOME)/include -e $(CSOURCE) > tmp.h
	mv tmp.h $@
