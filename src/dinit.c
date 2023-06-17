// dear emacs, please treat this as -*- C++ -*-

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Arith.h>

#include "pomp_internal.h"

static R_INLINE SEXP add_args
(
 SEXP args, SEXP Snames, SEXP Pnames, SEXP Cnames
 ) {

  SEXP var;
  int v;

  PROTECT(args);

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

// compute pdf of the initial state
static SEXP init_density
(
 SEXP func, SEXP X, SEXP t0, SEXP params, SEXP covar,
 SEXP log, SEXP args, SEXP gnsi
 ) {

  pompfunmode mode = undef;
  int give_log;
  int nvars, npars, nreps, ncovars;
  SEXP Snames, Pnames, Cnames;
  SEXP fn;
  SEXP F, cvec;
  double *cov;
  int *dim;

  dim = INTEGER(GET_DIM(X)); nvars = dim[0]; nreps = dim[1];
  dim = INTEGER(GET_DIM(params)); npars = dim[0];

  PROTECT(Snames = GET_ROWNAMES(GET_DIMNAMES(X)));
  PROTECT(Pnames = GET_ROWNAMES(GET_DIMNAMES(params)));
  PROTECT(Cnames = get_covariate_names(covar));

  PROTECT(F = ret_array(nreps));

  // set up the covariate table
  lookup_table_t covariate_table = make_covariate_table(covar,&ncovars);
  PROTECT(cvec = NEW_NUMERIC(ncovars));
  cov = REAL(cvec);

  PROTECT(fn = pomp_fun_handler(func,gnsi,&mode,Snames,Pnames,NA_STRING,Cnames));

  give_log = *(INTEGER(AS_INTEGER(log)));

  int nprotect = 6;

  switch (mode) {

  case Rfun: {

    SEXP ans;
    double *t = REAL(t0);
    double *ps;
    double *x = REAL(X);
    double *ft = REAL(F);
    int j;

    PROTECT(args = add_args(args,Snames,Pnames,Cnames)); nprotect++;

    // interpolate the covariates at time t1
    table_lookup(&covariate_table,*t,cov);

    for (j = 0, ps = REAL(params);
	 j < nreps;
	 j++, ft++, x += nvars, ps += npars) {

      PROTECT(ans = eval_call(fn,args,t,x,nvars,ps,npars,cov,ncovars));
      
      *ft = *(REAL(AS_NUMERIC(ans)));
      
      UNPROTECT(1);
      
      if (!give_log) *ft = exp(*ft);
      
    }


  }

    break;

  case native: case regNative: {

    int *sidx, *pidx, *cidx;
    double *t = REAL(t0);
    double *ps = REAL(params);
    double *x = REAL(X);
    double *ft = REAL(F);
    pomp_dinit *ff = NULL;
    int j;

    sidx = INTEGER(GET_SLOT(func,install("stateindex")));
    pidx = INTEGER(GET_SLOT(func,install("paramindex")));
    cidx = INTEGER(GET_SLOT(func,install("covarindex")));

    *((void **) (&ff)) = R_ExternalPtrAddr(fn);

    set_pomp_userdata(args);

    // interpolate the covariates
    table_lookup(&covariate_table,*t,cov);

    for (j = 0, ps = REAL(params);
	 j < nreps;
	 j++, ft++, x += nvars, ps += npars) {

      (*ff)(ft,x,ps,*t,sidx,pidx,cidx,cov);

      if (!give_log) *ft = exp(*ft);

    }

    unset_pomp_userdata();
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

  int *xdim, npars, nvars, nreps, nrepsx;
  SEXP F, fn, args, covar;

  PROTECT(t0=AS_NUMERIC(t0));

  PROTECT(x = as_matrix(x));
  xdim = INTEGER(GET_DIM(x));
  nvars = xdim[0]; nrepsx = xdim[1];

  PROTECT(params = as_matrix(params));
  xdim = INTEGER(GET_DIM(params));
  npars = xdim[0]; nreps = xdim[1];

  int nprotect = 3;

  if (nrepsx > nreps) {         // more states than parameters
    if (nrepsx % nreps != 0) {
      err("the larger number of replicates is not a multiple of smaller.");
    } else {
      SEXP copy;
      double *src, *tgt;
      int dims[2];
      int j, k;
      dims[0] = npars; dims[1] = nrepsx;
      PROTECT(copy = duplicate(params));
      PROTECT(params = makearray(2,dims));
      nprotect += 2;
      setrownames(params,GET_ROWNAMES(GET_DIMNAMES(copy)),2);
      src = REAL(copy);
      tgt = REAL(params);
      for (j = 0; j < nrepsx; j++) {
        for (k = 0; k < npars; k++, tgt++) {
          *tgt = src[k+npars*(j%nreps)];
        }
      }
    }
    nreps = nrepsx;
  } else if (nrepsx < nreps) {  // more parameters than states
    if (nreps % nrepsx != 0) {
      err("the larger number of replicates is not a multiple of smaller.");
    } else {
      SEXP copy;
      double *src, *tgt;
      int dims[2];
      int j, k;
      dims[0] = nvars; dims[1] = nreps;
      PROTECT(copy = duplicate(x));
      PROTECT(x = makearray(2,dims));
      nprotect += 2;
      setrownames(x,GET_ROWNAMES(GET_DIMNAMES(copy)),2);
      src = REAL(copy);
      tgt = REAL(x);
      for (j = 0; j < nreps; j++) {
	for (k = 0; k < nvars; k++, tgt++) {
	  *tgt = src[k+nvars*(j%nrepsx)];
	}
      }
    }
  }

  // extract the process function
  PROTECT(fn = GET_SLOT(object,install("dinit")));
  // extract other arguments
  PROTECT(args = VectorToPairList(GET_SLOT(object,install("userdata"))));
  PROTECT(covar = GET_SLOT(object,install("covar")));
  // evaluate the density
  PROTECT(F = init_density(fn,x,t0,params,covar,log,args,gnsi));

  nprotect += 4;

  UNPROTECT(nprotect);
  return F;
}
