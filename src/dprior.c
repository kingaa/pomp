// dear emacs, please treat this as -*- C++ -*-

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "pomp_internal.h"

static R_INLINE SEXP add_args (SEXP names, SEXP log, SEXP args)
{

  int nprotect = 0;
  SEXP var;
  int v;

  PROTECT(args = LCONS(AS_LOGICAL(log),args)); nprotect++;
  SET_TAG(args,install("log"));

  for (v = LENGTH(names)-1; v >= 0; v--) {
    PROTECT(var = NEW_NUMERIC(1)); nprotect++;
    PROTECT(args = LCONS(var,args)); nprotect++;
    SET_TAG(args,install(CHAR(STRING_ELT(names,v))));
  }

  UNPROTECT(nprotect);
  return args;

}

static R_INLINE SEXP eval_call (SEXP fn, SEXP args, double *p, int n)
{

  SEXP var = args, ans, ob;
  int v;

  for (v = 0; v < n; v++, p++, var=CDR(var))
    *(REAL(CAR(var))) = *p;

  PROTECT(ob = LCONS(fn,args));
  PROTECT(ans = eval(ob,CLOENV(fn)));

  UNPROTECT(2);
  return ans;

}

SEXP do_dprior (SEXP object, SEXP params, SEXP log, SEXP gnsi)
{
  int nprotect = 0;
  pompfunmode mode = undef;
  int npars, nreps;
  SEXP Pnames, pompfun, fn, args, F;
  int *dim;

  PROTECT(params = as_matrix(params)); nprotect++;
  dim = INTEGER(GET_DIM(params));
  npars = dim[0]; nreps = dim[1];

  PROTECT(Pnames = GET_ROWNAMES(GET_DIMNAMES(params))); nprotect++;

  // extract the user-defined function
  PROTECT(pompfun = GET_SLOT(object,install("dprior"))); nprotect++;
  PROTECT(fn = pomp_fun_handler(pompfun,gnsi,&mode,NA_STRING,Pnames,NA_STRING,NA_STRING)); nprotect++;

  // extract 'userdata' as pairlist
  PROTECT(args = VectorToPairList(GET_SLOT(object,install("userdata")))); nprotect++;

  // to store results
  PROTECT(F = NEW_NUMERIC(nreps)); nprotect++;

  switch (mode) {
  case Rfun: {
    SEXP ans;
    double *ps, *pt;
    int j;

    PROTECT(args = add_args(Pnames,log,args)); nprotect++;

    for (j = 0, ps = REAL(params), pt = REAL(F); j < nreps; j++, ps += npars, pt++) {

      PROTECT(ans = eval_call(fn,args,ps,npars));
      *pt = *(REAL(AS_NUMERIC(ans)));
      UNPROTECT(1);

    }
  }

    break;

  case native: case regNative: {
    int give_log, *pidx = 0;
    pomp_dprior *ff = NULL;
    double *ps, *pt;
    int j;

    // construct state, parameter, covariate, observable indices
    pidx = INTEGER(GET_SLOT(pompfun,install("paramindex")));

    // address of native routine
    *((void **) (&ff)) = R_ExternalPtrAddr(fn);

    give_log = *(INTEGER(AS_INTEGER(log)));

    R_CheckUserInterrupt();	// check for user interrupt

    set_pomp_userdata(args);

    // loop over replicates
    for (j = 0, pt = REAL(F), ps = REAL(params); j < nreps; j++, ps += npars, pt++)
      (*ff)(pt,ps,give_log,pidx);

    unset_pomp_userdata();
  }

    break;

  default: {
    int give_log, j;
    double *pt;

    give_log = *(INTEGER(AS_INTEGER(log)));

    // loop over replicates
    for (j = 0, pt = REAL(F); j < nreps; j++, pt++)
      *pt = (give_log) ? 0.0 : 1.0;

  }

  }

  UNPROTECT(nprotect);
  return F;
}
