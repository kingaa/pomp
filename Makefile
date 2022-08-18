include rules.mk

REVDEPS=spaero epimdr CollocInfer spatPomp DTAT

includes: inst/include/pomp.h inst/include/pomp_defines.h

headers: src/pomp_decls.h

src/pomp_decls.h: $(CSOURCE)
	cproto -I $(R_HOME)/include -e $(CSOURCE) > tmp.h
	mv tmp.h $@
