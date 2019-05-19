// dear emacs, please treat this as -*- C++ -*-

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <Rversion.h>
#include <stdarg.h>

#include "pomp_internal.h"

static R_INLINE SEXP name_index (SEXP provided, SEXP object, const char *slot, const char *humanreadable) {
  SEXP slotnames, index;
  PROTECT(slotnames = GET_SLOT(object,install(slot)));
  if (LENGTH(slotnames) > 0) {
    PROTECT(index = matchnames(provided,slotnames,humanreadable));
  } else {
    PROTECT(index = NEW_INTEGER(0));
  }
  UNPROTECT(2);
  return index;
}

// Returns either the R function or the address of the native routine.
// On return, mode indicates the mode of the 'pomp_fun'
// (i.e., R function, external function, or C snippet).
// If 'gnsi' is set to TRUE, we look up the native symbol information in the DLL,
// storing it in the 'address' slot.
// If 'gsni' is TRUE, and there are names in one or more of the S,P,O,C arguments, we look up the
// names in the corresponding 'pomp_fun' slots and storing the corresponding index
// inside the 'pomp_fun'.
SEXP pomp_fun_handler (SEXP pfun, SEXP gnsi, pompfunmode *mode,
  SEXP S, SEXP P, SEXP O, SEXP C)
{
  int nprotect = 0;
  SEXP f = R_NilValue;
  SEXP sidx, pidx, oidx, cidx;

  *mode = *(INTEGER(GET_SLOT(pfun,install("mode"))));

  switch (*mode) {

  case Rfun:			// R function

    PROTECT(f = GET_SLOT(pfun,install("R.fun"))); nprotect++;

    break;

  case native: case regNative:	// native code

    if (*(LOGICAL(gnsi))) {	// get native symbol information?

      SEXP nf, pack;
      PROTECT(nf = GET_SLOT(pfun,install("native.fun"))); nprotect++;
      PROTECT(pack = GET_SLOT(pfun,install("PACKAGE"))); nprotect++;
      if (LENGTH(pack) < 1) {
        PROTECT(pack = mkString("")); nprotect++;
      }

      switch (*mode) {
      case native: {
        SEXP nsi;
        PROTECT(nsi = eval(PROTECT(lang3(install("getNativeSymbolInfo"),nf,pack)),R_BaseEnv)); nprotect += 2;
        PROTECT(f = getListElement(nsi,"address")); nprotect++;
      }
        break;

      case regNative: {
        // Before version 3.4.0, R_MakeExternalPtrFn is not part of the R API.
        // Therefore, we must use some trickery to avoid the ISO C proscription of
        //     (void *) <-> (function *) conversion.
        const char *fname, *pkg;
        fname = (const char *) CHAR(STRING_ELT(nf,0));
        pkg = (const char *) CHAR(STRING_ELT(pack,0));
#if (R_VERSION < 197632) // before 3.4.0
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

      default: // #nocov
        break; // #nocov
      }

      SET_SLOT(pfun,install("address"),f);

      if (S != NA_STRING) {
        PROTECT(sidx = name_index(S,pfun,"statenames","state variables")); nprotect++;
        SET_SLOT(pfun,install("stateindex"),sidx);
      }

      if (P != NA_STRING) {
        PROTECT(pidx = name_index(P,pfun,"paramnames","parameters")); nprotect++;
        SET_SLOT(pfun,install("paramindex"),pidx);
      }

      if (O != NA_STRING) {
        PROTECT(oidx = name_index(O,pfun,"obsnames","observables")); nprotect++;
        SET_SLOT(pfun,install("obsindex"),oidx);
      }

      if (C != NA_STRING) {
        PROTECT(cidx = name_index(C,pfun,"covarnames","covariates")); nprotect++;
        SET_SLOT(pfun,install("covarindex"),cidx);
      }

    } else {			// native symbol info is stored

      PROTECT(f = GET_SLOT(pfun,install("address"))); nprotect++;

    }

    break;

  case undef: default: {

    PROTECT(f = R_NilValue); nprotect++;
    *mode = undef;

  }

  }

  UNPROTECT(nprotect);
  return f;
}

SEXP load_stack_incr (SEXP pack) {
  const char *pkg;
  void (*ff)(void);
  pkg = (const char *) CHAR(STRING_ELT(pack,0));
  ff = (void (*)(void)) R_GetCCallable(pkg,"__pomp_load_stack_incr");
  ff();
  return R_NilValue;
}

SEXP load_stack_decr (SEXP pack) {
  SEXP s;
  const char *pkg;
  void (*ff)(int *);
  PROTECT(s = NEW_INTEGER(1));
  pkg = (const char *) CHAR(STRING_ELT(pack,0));
  ff = (void (*)(int *)) R_GetCCallable(pkg,"__pomp_load_stack_decr");
  ff(INTEGER(s));
  if (*(INTEGER(s)) < 0) errorcall(R_NilValue,"impossible!");
  UNPROTECT(1);
  return s;
}

// SEXP concat (int nargs, ...) {
//   int nprotect = 0;
//   va_list ap;
//   SEXP f = R_NilValue;
//   int i;

//   va_start(ap,nargs);
//   for (i = 0; i < nargs; i++) {
//     PROTECT(f = LCONS(va_arg(ap,SEXP),f)); nprotect++;
//   }
//   va_end(ap);

//   PROTECT(f = eval(LCONS(install("c"),f),R_BaseEnv)); nprotect++;

//   UNPROTECT(nprotect);
//   return f;
// }

// 'pomp_fun' is provided for use by other packages
// It returns a list of two elements.
// The first is the the R function or the address of the native routine.
// The second is the 'mode'.

// SEXP pomp_fun (SEXP pfun, SEXP gnsi) {
//   pompfunmode md;
//   SEXP fun, mode, retval;
//   PROTECT(fun = pomp_fun_handler(pfun,gnsi,&md));
//   PROTECT(mode = NEW_INTEGER(1));
//   *(INTEGER(mode)) = (int) md;
//   PROTECT(retval = NEW_LIST(2));
//   SET_ELEMENT(retval,0,fun);
//   SET_ELEMENT(retval,0,mode);
//   return retval;
// }
