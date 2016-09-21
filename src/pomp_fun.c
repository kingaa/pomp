// dear emacs, please treat this as -*- C++ -*-

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <Rversion.h>

#include "pomp_internal.h"

// returns either the R function or the address of the native routine
// on return, use_native tells whether to use the native or the R function
SEXP pomp_fun_handler (SEXP pfun, SEXP gnsi, pompfunmode *mode)
{
  int nprotect = 0;
  SEXP f = R_NilValue;

  *mode = *(INTEGER(GET_SLOT(pfun,install("mode"))));

  switch (*mode) {

  case Rfun:			// R function

    PROTECT(f = GET_SLOT(pfun,install("R.fun"))); nprotect++;

    break;

  case native: case regNative:	// native code

    if (*(INTEGER(gnsi))) {	// get native symbol information?

      SEXP nf, pack;
      PROTECT(nf = GET_SLOT(pfun,install("native.fun"))); nprotect++;
      PROTECT(pack = GET_SLOT(pfun,install("PACKAGE"))); nprotect++;
      if (LENGTH(pack) < 1) {
	PROTECT(pack = mkString("")); nprotect++;
      }

      switch (*mode) {
      case native: 
	{
	  SEXP nsi;
	  PROTECT(nsi = eval(lang3(install("getNativeSymbolInfo"),nf,pack),R_BaseEnv)); nprotect++;
	  PROTECT(f = getListElement(nsi,"address")); nprotect++;
	}
	break;

      case regNative:
	{
	  // Before svn version 71180, R_MakeExternalPtrFn is not part of the R API.
	  // Therefore, we must use some trickery to avoid the ISO C proscription of
	  //     (void *) <-> (function *) conversion.
	  const char *fname, *pkg;
	  fname = (const char *) CHARACTER_DATA(STRING_ELT(nf,0));
	  pkg = (const char *) CHARACTER_DATA(STRING_ELT(pack,0));
	  int svn = R_SVN_REVISION;
#if (R_SVN_REVISION < 71180)
	  // This is cadged from 'R_MakeExternalPtrFn'.
	  union {void *p; DL_FUNC fn;} trick;
	  trick.fn = R_GetCCallable(pkg,fname);
	  PROTECT(f = R_MakeExternalPtr(trick.p,R_NilValue,R_NilValue)); nprotect++;
#else
	  DL_FUNC fn;
	  fn = R_GetCCallable(pkg,fname);
	  PROTECT(f = R_MakeExternalPtrFn(fn,R_NilValue,R_NilValue)); nprotect++;
#endif
	}
	break;
      
      case Rfun: case undef: default:
	break;			// # nocov
      }

      SET_SLOT(pfun,install("address"),f);

    } else {			// native symbol info is stored

      PROTECT(f = GET_SLOT(pfun,install("address"))); nprotect++;

    }

    *mode = native;

    break;

  case undef: default:
    {
      const char *purp = (const char *) CHARACTER_DATA(STRING_ELT(GET_SLOT(pfun,install("purpose")),0));

      errorcall(R_NilValue,"operation cannot be completed: %s has not been specified",purp);
    }

  }

  UNPROTECT(nprotect);
  return f;
}

SEXP load_stack_incr (SEXP pack) {
  const char *pkg;
  void (*ff)(void);
  pkg = (const char *) CHARACTER_DATA(STRING_ELT(pack,0));
  ff = (void (*)(void)) R_GetCCallable(pkg,"__pomp_load_stack_incr");
  ff();
  return R_NilValue;
}

SEXP load_stack_decr (SEXP pack) {
  SEXP s;
  const char *pkg;
  void (*ff)(int *);
  PROTECT(s = NEW_INTEGER(1));
  pkg = (const char *) CHARACTER_DATA(STRING_ELT(pack,0));
  ff = (void (*)(int *)) R_GetCCallable(pkg,"__pomp_load_stack_decr");
  ff(INTEGER(s));
  if (*(INTEGER(s)) < 0) errorcall(R_NilValue,"impossible!");
  UNPROTECT(1);
  return s;
}
