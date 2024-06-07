REPODIR = repo
INCLUDES = inst/include/pomp.h inst/include/pomp_defines.h
HEADERS = src/decls.h
SESSION_PKGS = datasets,utils,grDevices,graphics,stats,methods,tidyverse,pomp

include rules.mk

src/decls.h: $(CSOURCE)
	file=`mktemp tmpXXXXXXX.h` && \
	cproto -f2 -I $(R_HOME)/include -e $(CSOURCE) > $$file && \
	mv $$file $@
