// dear emacs, please treat this as -*- C++ -*-

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "pomp_internal.h"

static R_INLINE SEXP add_args (SEXP args, SEXP names)
{

  int nprotect = 0;
  SEXP var;
  int v;

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

  SEXP var = args, ans;
  int v;

  for (v = 0; v < n; v++, p++, var=CDR(var)) *(REAL(CAR(var))) = *p;

  PROTECT(ans = eval(LCONS(fn,args),CLOENV(fn)));

  UNPROTECT(1);
  return ans;

}

static R_INLINE SEXP ret_array (SEXP params)
{
  const char *dimnm[2] = {"variable", "rep"};
  SEXP P;

  PROTECT(P = duplicate(params));
  fixdimnames(P,dimnm,2);

  UNPROTECT(1);
  return P;

}

SEXP do_rprior (SEXP object, SEXP params, SEXP gnsi)
{
  int nprotect = 0;
  pompfunmode mode = undef;
  int npars, nreps;
  SEXP Pnames, pompfun, fn, args;
  SEXP P = R_NilValue;
  int *dim;

  PROTECT(params = as_matrix(params)); nprotect++;
  dim = INTEGER(GET_DIM(params));
  npars = dim[0]; nreps = dim[1];

  PROTECT(Pnames = GET_ROWNAMES(GET_DIMNAMES(params))); nprotect++;

  // extract the user-defined function
  PROTECT(pompfun = GET_SLOT(object,install("rprior"))); nprotect++;
  PROTECT(fn = pomp_fun_handler(pompfun,gnsi,&mode)); nprotect++;

  // extract 'userdata' as pairlist
  PROTECT(args = VectorToPairList(GET_SLOT(object,install("userdata")))); nprotect++;
  PROTECT(P = ret_array(params)); nprotect++;

  switch (mode) {

  case Rfun: {

    SEXP ans, nm;
    double *pa, *ps = REAL(params), *pt = REAL(P);
    int *posn = 0;
    int i, j;

    // set up the function call
    PROTECT(args = add_args(args,Pnames)); nprotect++;

    for (j = 0; j < nreps; j++, ps += npars, pt += npars) {

      if (j == 0) {

        PROTECT(ans = eval_call(fn,args,ps,npars)); nprotect++;

        PROTECT(nm = GET_NAMES(ans)); nprotect++;
        if (isNull(nm))
          errorcall(R_NilValue,"'rprior' must return a named numeric vector.");
        posn = INTEGER(PROTECT(matchnames(Pnames,nm,"parameters"))); nprotect++;

        pa = REAL(AS_NUMERIC(ans));

        for (i = 0; i < LENGTH(ans); i++) pt[posn[i]] = pa[i];

      } else {

        PROTECT(ans = eval_call(fn,args,ps,npars));
        pa = REAL(AS_NUMERIC(ans));
        for (i = 0; i < LENGTH(ans); i++) pt[posn[i]] = pa[i];
        UNPROTECT(1);

      }
    }

  }

    break;

  case native: {

    double *ps;
    int *pidx = 0;
    pomp_rprior *ff = NULL;
    int j;

    // construct parameter indices
    pidx = INTEGER(PROTECT(name_index(Pnames,pompfun,"paramnames","parameters"))); nprotect++;

    // address of native routine
    *((void **) (&ff)) = R_ExternalPtrAddr(fn);

    R_CheckUserInterrupt();	// check for user interrupt

    set_pomp_userdata(args);
    GetRNGstate();

    // loop over replicates
    for (j = 0, ps = REAL(P); j < nreps; j++, ps += npars)
      (*ff)(ps,pidx);

    PutRNGstate();
    unset_pomp_userdata();
  }

    break;

  default: // just duplicate

    break;

  }

  UNPROTECT(nprotect);
  return P;
}
