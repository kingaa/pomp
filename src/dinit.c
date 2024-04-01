// dear emacs, please treat this as -*- C++ -*-

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

#include "internal.h"

static R_INLINE SEXP add_args
(
 SEXP args, SEXP Snames, SEXP Pnames, SEXP Cnames
 ) {

  SEXP var;
  int v;

  PROTECT(args = VectorToPairList(args));

  // Covariates
  for (v = LENGTH(Cnames)-1; v >= 0; v--) {
    var = NEW_NUMERIC(1);
    args = LCONS(var,args);
    UNPROTECT(1);
    PROTECT(args);
    SET_TAG(args,installChar(STRING_ELT(Cnames,v)));
  }

  // Parameters
  for (v = LENGTH(Pnames)-1; v >= 0; v--) {
    var = NEW_NUMERIC(1);
    args = LCONS(var,args);
    UNPROTECT(1);
    PROTECT(args);
    SET_TAG(args,installChar(STRING_ELT(Pnames,v)));
  }

  // Latent state variables
  for (v = LENGTH(Snames)-1; v >= 0; v--) {
    var = NEW_NUMERIC(1);
    args = LCONS(var,args);
    UNPROTECT(1);
    PROTECT(args);
    SET_TAG(args,installChar(STRING_ELT(Snames,v)));
  }

  // Time
  var = NEW_NUMERIC(1);
  args = LCONS(var,args);
  UNPROTECT(1);
  PROTECT(args);
  SET_TAG(args,install("t0"));

  UNPROTECT(1);
  return args;

}

static R_INLINE SEXP eval_call
(
 SEXP fn, SEXP args,
 double *t0,
 double *x, int nvar,
 double *p, int npar,
 double *c, int ncov
 ) {

  SEXP var = args, ans, ob;
  int v;

  *(REAL(CAR(var))) = *t0; var = CDR(var);
  for (v = 0; v < nvar; v++, x++) {
    *(REAL(CAR(var))) = *x; var = CDR(var);
  }
  for (v = 0; v < npar; v++, p++, var=CDR(var)) *(REAL(CAR(var))) = *p;
  for (v = 0; v < ncov; v++, c++, var=CDR(var)) *(REAL(CAR(var))) = *c;

  PROTECT(ob = LCONS(fn,args));
  PROTECT(ans = eval(ob,CLOENV(fn)));

  UNPROTECT(2);
  return ans;

}

static R_INLINE SEXP ret_array (int nreps)
{
  int dim[1] = {nreps};
  const char *dimnm[1] = {".id"};
  SEXP F;
  PROTECT(F = makearray(1,dim));
  fixdimnames(F,dimnm,1);
  UNPROTECT(1);
  return F;
}

static SEXP init_density
(
 SEXP func, SEXP X, SEXP t0, SEXP params, SEXP covar,
 SEXP log, SEXP args, SEXP gnsi
 ) {

  pompfunmode mode = undef;
  int give_log;
  int nvars, npars, nrepsx, nrepsp, nreps, ncovars;
  SEXP Snames, Pnames, Cnames;
  SEXP fn, ans;
  SEXP F, cvec;
  double *cov;
  int *dim;

  dim = INTEGER(GET_DIM(X)); nvars = dim[0]; nrepsx = dim[1];
  dim = INTEGER(GET_DIM(params)); npars = dim[0]; nrepsp = dim[1];

  give_log = *(INTEGER(AS_INTEGER(log)));

  // handle case with different numbers of states and parameters
  if (nrepsx != nrepsp && nrepsx % nrepsp != 0 && nrepsp % nrepsx != 0) {
    err("the larger number of replicates is not a multiple of smaller.");
  } else {
    nreps = (nrepsx > nrepsp) ? nrepsx : nrepsp;
  }

  PROTECT(Snames = GET_ROWNAMES(GET_DIMNAMES(X)));
  PROTECT(Pnames = GET_ROWNAMES(GET_DIMNAMES(params)));
  PROTECT(Cnames = get_covariate_names(covar));
  PROTECT(F = ret_array(nreps));

  // set up the covariate table
  lookup_table_t covariate_table = make_covariate_table(covar,&ncovars);
  PROTECT(cvec = NEW_NUMERIC(ncovars));
  cov = REAL(cvec);

  // process the pomp_fun
  PROTECT(fn = pomp_fun_handler(func,gnsi,&mode,Snames,Pnames,NA_STRING,Cnames));

  int nprotect = 6;

  switch (mode) {

  case Rfun: {

    double *t = REAL(t0);
    double *ft = REAL(F);

    PROTECT(args = add_args(args,Snames,Pnames,Cnames)); nprotect++;

    // interpolate the covariates at time t1
    table_lookup(&covariate_table,*t,cov);

    for (int j = 0; j < nreps; j++, ft++) {

      double *xs = REAL(X)+nvars*(j%nrepsx);
      double *ps = REAL(params)+npars*(j%nrepsp);

      PROTECT(ans = eval_call(fn,args,t,xs,nvars,ps,npars,cov,ncovars));

      *ft = *REAL(AS_NUMERIC(ans));

      UNPROTECT(1);

      if (!give_log) *ft = exp(*ft);

    }


  }

    break;

  case native: case regNative: {

    int *sidx, *pidx, *cidx;
    double *t = REAL(t0);
    double *ft = REAL(F);
    pomp_dinit *ff = NULL;

    sidx = INTEGER(GET_SLOT(func,install("stateindex")));
    pidx = INTEGER(GET_SLOT(func,install("paramindex")));
    cidx = INTEGER(GET_SLOT(func,install("covarindex")));

    *((void **) (&ff)) = R_ExternalPtrAddr(fn);

    // interpolate the covariates
    table_lookup(&covariate_table,*t,cov);

    for (int j = 0; j < nreps; j++, ft++) {

      double *xs = REAL(X)+nvars*(j%nrepsx);
      double *ps = REAL(params)+npars*(j%nrepsp);

      (*ff)(ft,xs,ps,*t,sidx,pidx,cidx,cov);

      if (!give_log) *ft = exp(*ft);

    }

  }

    break;

  default: {
    double *ft = REAL(F);
    int j;

    for (j = 0; j < nreps; j++, ft++) { // loop over replicates
      *ft = R_NaReal;
    }

    warn("'dinit' unspecified: likelihood undefined.");

  }

  }

  UNPROTECT(nprotect);
  return F;
}

SEXP do_dinit
(
 SEXP object, SEXP t0, SEXP x, SEXP params, SEXP log, SEXP gnsi
 ) {
  SEXP F, fn, args, covar;
  PROTECT(t0=AS_NUMERIC(t0));
  PROTECT(x = as_matrix(x));
  PROTECT(params = as_matrix(params));
  // extract the process function
  PROTECT(fn = GET_SLOT(object,install("dinit")));
  // extract other arguments
  PROTECT(args = GET_SLOT(object,install("userdata")));
  PROTECT(covar = GET_SLOT(object,install("covar")));
  // evaluate the density
  PROTECT(F = init_density(fn,x,t0,params,covar,log,args,gnsi));
  UNPROTECT(7);
  return F;
}
