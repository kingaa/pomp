include rules.mk

REVDEPS=spaero epimdr CollocInfer spatPomp DTAT
TOKEN=cbbc302d-fff3-4530-8474-0f3f48db6776

includes: inst/include/pomp.h inst/include/pomp_defines.h

headers: src/pomp_decls.h

src/pomp_decls.h: $(CSOURCE)
	cproto -I $(R_HOME)/include -e $(CSOURCE) > tmp.h
	mv tmp.h $@
