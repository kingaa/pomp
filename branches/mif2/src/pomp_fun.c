// dear emacs, please treat this as -*- C++ -*-

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "pomp_internal.h"

// returns either the R function or the address of the native routine
// on return, use_native tells whether to use the native or the R function
SEXP pomp_fun_handler (SEXP pfun, SEXP gnsi, int *mode) 
{
  int nprotect = 0;
  SEXP nf, pack, nsi, f = R_NilValue;

  *mode = *(INTEGER(GET_SLOT(pfun,install("mode"))))-1;

  switch (*mode) {
  case 0:			// R function
    PROTECT(f = GET_SLOT(pfun,install("R.fun"))); nprotect++;
    break;
  case 1:			// native code
    if (*(INTEGER(gnsi))) {	// get native symbol information
      PROTECT(nf = GET_SLOT(pfun,install("native.fun"))); nprotect++;
      PROTECT(pack = GET_SLOT(pfun,install("PACKAGE"))); nprotect++;
      if (LENGTH(pack) < 1) {
	PROTECT(pack = NEW_CHARACTER(1)); nprotect++;
	SET_STRING_ELT(pack,0,mkChar(""));
      }
      PROTECT(nsi = eval(lang3(install("getNativeSymbolInfo"),nf,pack),R_BaseEnv)); nprotect++;
      PROTECT(f = getListElement(nsi,"address")); nprotect++;
      SET_SLOT(pfun,install("address"),f);
    } else {
      PROTECT(f = GET_SLOT(pfun,install("address"))); nprotect++;
    }
    break;
  default:
    error("operation cannot be completed: some needed function has not been specified");
    break;
  }

  UNPROTECT(nprotect);
  return f;
}
