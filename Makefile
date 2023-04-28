REPODIR = repo
INCLUDES = inst/include/pomp.h inst/include/pomp_defines.h
HEADERS = src/pomp_decls.h

include rules.mk

src/pomp_decls.h: $(CSOURCE)
	file=`mktemp tmpXXXXXXX.h` && \
	cproto -I $(R_HOME)/include -e $(CSOURCE) > $$file && \
	mv $$file $@
