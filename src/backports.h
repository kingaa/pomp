// -*- C++ -*-

#ifndef _POMP_BACKPORTS_H_
#define _POMP_BACKPORTS_H_

#include <Rversion.h>

// Some backports needed for older versions of R.
// Taken from
// https://cran.r-project.org/doc/manuals/r-devel/R-exts.html#Some-backports
// These are not necessarily all used in the package code.

#if R_VERSION < R_Version(4, 5, 0)

# define isDataFrame(x) Rf_isFrame(x)
# define R_ClosureFormals(x) FORMALS(x)
# define R_ClosureEnv(x) CLOENV(x)
# define R_ParentEnv(x) ENCLOS(x)

#endif

#endif
