// dear emacs, please treat this as -*- C++ -*-

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "pomp_internal.h"

// returns either the R function or the address of the native routine
// on return, use_native tells whether to use the native or the R function
SEXP pomp_fun_handler (SEXP pfun, int *mode) 
{
  int nprotect = 0;
  SEXP nf, pack, nsi, f = R_NilValue;
  int mod;

  mod = INTEGER(GET_SLOT(pfun,install("mode")))[0];
  *mode = mod-1;

  switch (*mode) {
  case 0:			// R function
    PROTECT(f = GET_SLOT(pfun,install("R.fun"))); nprotect++;
    break;
  case 1:			// native code
    PROTECT(nf = GET_SLOT(pfun,install("native.fun"))); nprotect++;
    PROTECT(pack = GET_SLOT(pfun,install("PACKAGE"))); nprotect++;
    if (LENGTH(pack) < 1) {
      PROTECT(pack = NEW_CHARACTER(1)); nprotect++;
      SET_STRING_ELT(pack,0,mkChar(""));
    }
    PROTECT(nsi = eval(lang3(install("getNativeSymbolInfo"),nf,pack),R_BaseEnv)); nprotect++;
    PROTECT(f = getListElement(nsi,"address")); nprotect++;
    break;
  default:
    error("'pomp_fun_handler': invalid 'mode' value");
    break;
  }

  UNPROTECT(nprotect);
  return f;
}

// returns either the R function or the address of the native routine
// returns a list of length 2
SEXP get_pomp_fun (SEXP pfun) {
  int nprotect = 0;
  SEXP ans, mode, f = R_NilValue;
  PROTECT(mode = NEW_INTEGER(1)); nprotect++;
  PROTECT(f = pomp_fun_handler(pfun,INTEGER(mode))); nprotect++;
  PROTECT(ans = NEW_LIST(2)); nprotect++;
  SET_ELEMENT(ans,0,f);	    // R function or pointer to native routine
  SET_ELEMENT(ans,1,mode);   // =0 if R function; =1 if native
  UNPROTECT(nprotect);
  return ans;
}

SEXP unpack_pomp_fun (SEXP pfunlist, int *mode) {
  int nprotect = 0;
  SEXP f;
  PROTECT(f = VECTOR_ELT(pfunlist,0)); nprotect++;
  *mode = INTEGER(VECTOR_ELT(pfunlist,1))[0];
  UNPROTECT(nprotect);
  return f;
}
