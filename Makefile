REPODIR = repo
INCLUDES = inst/include/pomp.h inst/include/pomp_defines.h
HEADERS = src/decls.h

include rules.mk

src/decls.h: $(CSOURCE)
	file=`mktemp tmpXXXXXXX.h` && \
	cproto -I $(R_HOME)/include -e $(CSOURCE) > $$file && \
	mv $$file $@
