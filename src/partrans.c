// dear emacs, please treat this as -*- C++ -*-

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <string.h>

#include "pomp_internal.h"

typedef enum {to = 1, from = -1} direction_t;

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

  SEXP var = args, ans, ob;
  int v;

  for (v = 0; v < n; v++, p++, var=CDR(var)) *(REAL(CAR(var))) = *p;

  PROTECT(ob = LCONS(fn,args));
  PROTECT(ans = eval(ob,CLOENV(fn)));

  UNPROTECT(2);
  return ans;

}

SEXP do_partrans (SEXP object, SEXP params, SEXP dir, SEXP gnsi)
{
  int nprotect = 0;
  SEXP Pnames, tparams, pompfun, fn, args, ob;
  pompfunmode mode = undef;
  direction_t direc;
  int qvec, npars, nreps;
  int *dim;

  qvec = isNull(GET_DIM(params)); // is 'params' a vector?

  PROTECT(tparams = duplicate(params)); nprotect++;

  // coerce 'params' to matrix
  PROTECT(tparams = as_matrix(tparams)); nprotect++;
  dim = INTEGER(GET_DIM(tparams));
  npars = dim[0]; nreps = dim[1];

  PROTECT(Pnames = GET_ROWNAMES(GET_DIMNAMES(tparams))); nprotect++;

  // determine direction of transformation and extract corresponding pomp_fun
  direc = (direction_t) *(INTEGER(dir));
  PROTECT(ob = GET_SLOT(object,install("partrans"))); nprotect++;
  switch (direc) {
  case from: default:	// from estimation scale
    PROTECT(pompfun = GET_SLOT(ob,install("from"))); nprotect++;
    break;
  case to:			// to estimation scale
    PROTECT(pompfun = GET_SLOT(ob,install("to"))); nprotect++;
    break;
  }

  PROTECT(fn = pomp_fun_handler(pompfun,gnsi,&mode,NA_STRING,Pnames,NA_STRING,NA_STRING)); nprotect++;

  // extract 'userdata' as pairlist
  PROTECT(args = VectorToPairList(GET_SLOT(object,install("userdata")))); nprotect++;

  switch (mode) {

  case Rfun: {

    SEXP ans, nm;
    double *pa, *ps = REAL(tparams);
    int *posn;
    int i, j;

    PROTECT(args = add_args(args,Pnames)); nprotect++;

    PROTECT(ans = eval_call(fn,args,ps,npars)); nprotect++;

    PROTECT(nm = GET_NAMES(ans)); nprotect++;
    if (invalid_names(nm))
      errorcall(R_NilValue,"user transformation functions must return named numeric vectors.");
    posn = INTEGER(PROTECT(matchnames(Pnames,nm,"parameters"))); nprotect++;

    pa = REAL(AS_NUMERIC(ans));

    for (i = 0; i < LENGTH(ans); i++) ps[posn[i]] = pa[i];

    for (j = 1, ps += npars; j < nreps; j++, ps += npars) {

      PROTECT(ans = eval_call(fn,args,ps,npars));
      pa = REAL(AS_NUMERIC(ans));
      for (i = 0; i < LENGTH(ans); i++) ps[posn[i]] = pa[i];
      UNPROTECT(1);

    }

  }

    break;

  case native: case regNative: {

    pomp_transform_fn *ff;
    double *ps, *pt;
    int *idx;
    int j;

    *((void **) (&ff)) = R_ExternalPtrAddr(fn);

    R_CheckUserInterrupt();

    idx = INTEGER(GET_SLOT(pompfun,install("paramindex")));

    set_pomp_userdata(args);

    for (j = 0, ps = REAL(params), pt = REAL(tparams); j < nreps; j++, ps += npars, pt += npars)
      (*ff)(pt,ps,idx);

    unset_pomp_userdata();

  }

    break;

  default:  // #nocov

    break;  // #nocov

  }

  if (qvec) {
    PROTECT(Pnames = GET_ROWNAMES(GET_DIMNAMES(tparams))); nprotect++;
    SET_DIM(tparams,R_NilValue);
    SET_NAMES(tparams,Pnames);
  }

  UNPROTECT(nprotect);
  return tparams;

}
