// dear emacs, please treat this as -*- C++ -*-

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <Rversion.h>
#include <stdarg.h>

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
        PROTECT(nsi = eval(PROTECT(lang3(install("getNativeSymbolInfo"),nf,pack)),R_BaseEnv)); nprotect += 2;
        PROTECT(f = getListElement(nsi,"address")); nprotect++;
      }
        break;

      case regNative:
      {
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

SEXP pomp_fun_args (SEXP args, SEXP Onames, SEXP Snames, SEXP Pnames, SEXP Cnames)
{
  int nprotect = 0;
  SEXP var;
  int v;

  // we construct the call from end to beginning
  // covariates, parameter, states, observables, then time

  // Covariates
  for (v = LENGTH(Cnames)-1; v >= 0; v--) {
    PROTECT(var = NEW_NUMERIC(1)); nprotect++;
    PROTECT(args = LCONS(var,args)); nprotect++;
    SET_TAG(args,install(CHAR(STRING_ELT(Cnames,v))));
  }

  // Parameters
  for (v = LENGTH(Pnames)-1; v >= 0; v--) {
    PROTECT(var = NEW_NUMERIC(1)); nprotect++;
    PROTECT(args = LCONS(var,args)); nprotect++;
    SET_TAG(args,install(CHAR(STRING_ELT(Pnames,v))));
  }

  // Latent state variables
  for (v = LENGTH(Snames)-1; v >= 0; v--) {
    PROTECT(var = NEW_NUMERIC(1)); nprotect++;
    PROTECT(args = LCONS(var,args)); nprotect++;
    SET_TAG(args,install(CHAR(STRING_ELT(Snames,v))));
  }

  // Observables
  for (v = LENGTH(Onames)-1; v >= 0; v--) {
    PROTECT(var = NEW_NUMERIC(1)); nprotect++;
    PROTECT(args = LCONS(var,args)); nprotect++;
    SET_TAG(args,install(CHAR(STRING_ELT(Onames,v))));
  }

  // Time
  PROTECT(var = NEW_NUMERIC(1)); nprotect++;
  PROTECT(args = LCONS(var,args)); nprotect++;
  SET_TAG(args,install("t"));

  UNPROTECT(nprotect);
  return args;

}

SEXP eval_pomp_fun_R_call (
    SEXP fn, SEXP args,
    double *t,
    double *y, int nobs,
    double *x, int nvar,
    double *p, int npar,
    double *c, int ncov)
{

  SEXP var = args, ans;
  int v;

  *(REAL(CAR(var))) = *t; var = CDR(var);
  for (v = 0; v < nobs; v++, y++, var=CDR(var)) *(REAL(CAR(var))) = *y;
  for (v = 0; v < nvar; v++, x++, var=CDR(var)) *(REAL(CAR(var))) = *x;
  for (v = 0; v < npar; v++, p++, var=CDR(var)) *(REAL(CAR(var))) = *p;
  for (v = 0; v < ncov; v++, c++, var=CDR(var)) *(REAL(CAR(var))) = *c;

  PROTECT(ans = eval(LCONS(fn,args),CLOENV(fn)));

  UNPROTECT(1);
  return ans;

}

SEXP pomp_fun_indices (SEXP fn, SEXP Onames, SEXP Snames, SEXP Pnames, SEXP Cnames)
{
  SEXP O = R_NilValue, S = R_NilValue, P = R_NilValue, C = R_NilValue;
  SEXP ans;

  // construct state, parameter, covariate, observable indices
  PROTECT(O = name_index(Onames,fn,"obsnames","observables"));
  PROTECT(S = name_index(Snames,fn,"statenames","state variables"));
  PROTECT(P = name_index(Pnames,fn,"paramnames","parameters"));
  PROTECT(C = name_index(Cnames,fn,"covarnames","covariates"));

  PROTECT(ans = NEW_LIST(4));
  SET_ELEMENT(ans,0,O);
  SET_ELEMENT(ans,1,S);
  SET_ELEMENT(ans,2,P);
  SET_ELEMENT(ans,3,C);

  UNPROTECT(5);
  return ans;

}

SEXP concat (int nargs, ...) {
  int nprotect = 0;
  va_list ap;
  SEXP f = R_NilValue;
  int i;

  va_start(ap,nargs);
  for (i = 0; i < nargs; i++) {
    PROTECT(f = LCONS(va_arg(ap,SEXP),f)); nprotect++;
  }
  va_end(ap);

  PROTECT(f = eval(LCONS(install("c"),f),R_BaseEnv)); nprotect++;

  UNPROTECT(nprotect);
  return f;
}
