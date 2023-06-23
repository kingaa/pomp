// dear emacs, please treat this as -*- C++ -*-

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

#include "internal.h"

typedef enum {to = 1, from = -1} direction_t;

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

SEXP do_partrans (SEXP object, SEXP params, SEXP dir, SEXP gnsi)
{

  SEXP Pnames, tparams, pompfun, fn, args, ob;
  pompfunmode mode = undef;
  direction_t direc;
  int qvec, npars, nreps;
  int *dim;

  qvec = isNull(GET_DIM(params)); // is 'params' a vector?

  PROTECT(tparams = duplicate(params));

  // coerce 'params' to matrix
  PROTECT(tparams = as_matrix(tparams));
  dim = INTEGER(GET_DIM(tparams));
  npars = dim[0]; nreps = dim[1];

  PROTECT(Pnames = GET_ROWNAMES(GET_DIMNAMES(tparams)));

  // determine direction of transformation and extract corresponding pomp_fun
  direc = (direction_t) *(INTEGER(dir));
  PROTECT(ob = GET_SLOT(object,install("partrans")));
  switch (direc) {
  case from: default:	// from estimation scale
    PROTECT(pompfun = GET_SLOT(ob,install("from")));
    break;
  case to:			// to estimation scale
    PROTECT(pompfun = GET_SLOT(ob,install("to")));
    break;
  }

  PROTECT(fn = pomp_fun_handler(pompfun,gnsi,&mode,NA_STRING,Pnames,NA_STRING,NA_STRING));

  // extract 'userdata' as pairlist
  PROTECT(args = VectorToPairList(GET_SLOT(object,install("userdata"))));

  int nprotect = 7;

  switch (mode) {

  case Rfun: {

    SEXP ans, nm;
    double *pa, *ps = REAL(tparams);
    int *posn;
    int i, j;

    PROTECT(args = add_args(args,Pnames));
    PROTECT(ans = eval_call(fn,args,ps,npars));

    PROTECT(nm = GET_NAMES(ans));
    if (invalid_names(nm))
      err("user transformation functions must return named numeric vectors.");
    posn = INTEGER(PROTECT(matchnames(Pnames,nm,"parameters")));

    nprotect += 4;

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

    pomp_transform *ff;
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
