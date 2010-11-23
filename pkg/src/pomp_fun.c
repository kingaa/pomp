// dear emacs, please treat this as -*- C++ -*-

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "pomp_internal.h"

// returns either the R function or the address of the native routine
// on return, use_native tells whether to use the native or the R function
SEXP pomp_fun_handler (SEXP pfun, int *use_native) 
{
  int nprotect = 0;
  SEXP nf, pack, nsi, f = R_NilValue;
  int use;

  use = INTEGER(GET_SLOT(pfun,install("use")))[0];

  switch (use) {
  case 1:			// use R function
    PROTECT(f = GET_SLOT(pfun,install("R.fun"))); nprotect++;
    *use_native = 0;
    break;
  case 2:			// use native routine
    PROTECT(nf = GET_SLOT(pfun,install("native.fun"))); nprotect++;
    PROTECT(pack = GET_SLOT(pfun,install("PACKAGE"))); nprotect++;
    if (LENGTH(pack) < 1) {
      PROTECT(pack = NEW_CHARACTER(1)); nprotect++;
      SET_STRING_ELT(pack,0,mkChar(""));
    }
    PROTECT(nsi = eval(lang3(install("getNativeSymbolInfo"),nf,pack),R_BaseEnv)); nprotect++;
    *use_native = 1;
    PROTECT(f = VECTOR_ELT(nsi,1)); nprotect++; // return a pointer to the loaded routine
    break;
  default:
    error("'pomp_fun_handler': invalid 'use' value");
    break;
  }

  UNPROTECT(nprotect);
  return f;
}

// returns either the R function or the address of the native routine
// returns a list of length 2
SEXP get_pomp_fun (SEXP pfun) {
  int nprotect = 0;
  SEXP ans, use, f = R_NilValue;
  PROTECT(use = NEW_INTEGER(1)); nprotect++;
  PROTECT(f = pomp_fun_handler(pfun,INTEGER(use))); nprotect++;
  PROTECT(ans = NEW_LIST(2)); nprotect++;
  SET_ELEMENT(ans,0,f);	    // R function or pointer to native routine
  SET_ELEMENT(ans,1,use);   // =0 if R function; =1 if native
  UNPROTECT(nprotect);
  return ans;
}
