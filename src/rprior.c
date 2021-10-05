// dear emacs, please treat this as -*- C++ -*-

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "pomp_internal.h"

static R_INLINE SEXP add_args (SEXP args, SEXP names)
{

  SEXP var;
  int v;

  PROTECT(args);

  for (v = LENGTH(names)-1; v >= 0; v--) {
    var = NEW_NUMERIC(1);
    args = LCONS(var,args);
    UNPROTECT(1);
    PROTECT(args);
    SET_TAG(args,installChar(STRING_ELT(names,v)));
  }

  UNPROTECT(1);
  return args;

}

static R_INLINE SEXP eval_call (SEXP fn, SEXP args, double *p, int n)
{

  SEXP var = args, ans, ob;
  int v;

  for (v = 0; v < n; v++, p++, var=CDR(var)) *(REAL(CAR(var))) = *p;

  PROTECT(ob = LCONS(fn,args));
  PROTECT(ans = eval(ob,CLOENV(fn)));

  UNPROTECT(2);
  return ans;

}

static R_INLINE SEXP ret_array (SEXP params)
{
  const char *dimnm[2] = {"variable", ".id"};
  SEXP P;

  PROTECT(P = duplicate(params));
  fixdimnames(P,dimnm,2);

  UNPROTECT(1);
  return P;

}

SEXP do_rprior (SEXP object, SEXP params, SEXP gnsi)
{

  pompfunmode mode = undef;
  int npars, nreps;
  SEXP Pnames, pompfun, fn, args;
  SEXP P = R_NilValue;
  int *dim;
  
  PROTECT(params = as_matrix(params));
  dim = INTEGER(GET_DIM(params));
  npars = dim[0]; nreps = dim[1];

  PROTECT(Pnames = GET_ROWNAMES(GET_DIMNAMES(params)));

  // extract the user-defined function
  PROTECT(pompfun = GET_SLOT(object,install("rprior")));
  PROTECT(fn = pomp_fun_handler(pompfun,gnsi,&mode,NA_STRING,Pnames,NA_STRING,NA_STRING));

  // extract 'userdata' as pairlist
  PROTECT(args = VectorToPairList(GET_SLOT(object,install("userdata"))));
  PROTECT(P = ret_array(params));

  int nprotect = 6;
  int first = 1;
  
  switch (mode) {
    
  case Rfun: {

    SEXP ans, nm;
    double *pa, *ps = REAL(params), *pt = REAL(P);
    int *posn = NULL;
    int i, j;

    // set up the function call
    PROTECT(args = add_args(args,Pnames)); nprotect++;

    for (j = 0; j < nreps; j++, ps += npars, pt += npars) {

      if (first) {

        PROTECT(ans = eval_call(fn,args,ps,npars));

        PROTECT(nm = GET_NAMES(ans));
        if (invalid_names(nm))
          err("'rprior' must return a named numeric vector.");
        posn = INTEGER(PROTECT(matchnames(Pnames,nm,"parameters")));
	
	nprotect += 3;

        pa = REAL(AS_NUMERIC(ans));
        for (i = 0; i < LENGTH(ans); i++) pt[posn[i]] = pa[i];

	first = 0;

      } else {

        PROTECT(ans = eval_call(fn,args,ps,npars));

        pa = REAL(AS_NUMERIC(ans));
        for (i = 0; i < LENGTH(ans); i++) pt[posn[i]] = pa[i];

        UNPROTECT(1);

      }
    }
  }

    break;

  case native: case regNative: {

    double *ps;
    int *pidx = 0;
    pomp_rprior *ff = NULL;
    int j;

    // extract parameter index
    pidx = INTEGER(GET_SLOT(pompfun,install("paramindex")));

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

    warn("'rprior' unspecified: duplicating parameters.");

  }

  UNPROTECT(nprotect);
  return P;
}
