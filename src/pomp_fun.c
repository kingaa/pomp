// dear emacs, please treat this as -*- C++ -*-

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>

#include "pomp_internal.h"

// returns either the R function or the address of the native routine
// on return, use_native tells whether to use the native or the R function

SEXP pomp_fun_handler (SEXP pfun, int *use_native) 
{
  int nprotect = 0;
  SEXP nsi, f = R_NilValue;
  int use = INTEGER(GET_SLOT(pfun,install("use")))[0];

  switch (use) {
  case 1:			// use R function
    PROTECT(f = GET_SLOT(pfun,install("R.fun"))); nprotect++;
    *use_native = 0;
    break;
  case 2:			// use native routine
    PROTECT(nsi = eval(
		       lang3(
			     install("getNativeSymbolInfo"),
			     GET_SLOT(pfun,install("native.fun")),
			     GET_SLOT(pfun,install("PACKAGE"))
			     ),
		       R_GlobalEnv
		       )
	    ); nprotect++;
    *use_native = 1;
    PROTECT(f = VECTOR_ELT(nsi,1)); nprotect++;
    break;
  default:
    UNPROTECT(nprotect);
    error("pomp_fun_handler: invalid value of 'use'");
    break;
  }

  UNPROTECT(nprotect);
  return f;
}
