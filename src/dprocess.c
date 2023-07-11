// dear emacs, please treat this as -*- C++ -*-

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

#include "internal.h"

static R_INLINE SEXP paste0 (SEXP a, SEXP b) {
  return eval(lang3(install("paste0"),a,b),R_BaseEnv);
}

static R_INLINE SEXP add_args (SEXP args, SEXP Snames, SEXP Pnames, SEXP Cnames)
{

  SEXP S1names, S2names;
  SEXP var;
  int v;

  PROTECT(S1names = paste0(Snames,mkString("_1")));
  PROTECT(S2names = paste0(Snames,mkString("_2")));

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
    SET_TAG(args,installChar(STRING_ELT(S2names,v)));

    var = NEW_NUMERIC(1);
    args = LCONS(var,args);
    UNPROTECT(1);
    PROTECT(args);
    SET_TAG(args,installChar(STRING_ELT(S1names,v)));

  }

  // Time
  var = NEW_NUMERIC(1);
  args = LCONS(var,args);
  UNPROTECT(1);
  PROTECT(args);
  SET_TAG(args,install("t_2"));

  var = NEW_NUMERIC(1);
  args = LCONS(var,args);
  UNPROTECT(1);
  PROTECT(args);
  SET_TAG(args,install("t_1"));

  UNPROTECT(3);
  return args;

}

static R_INLINE SEXP eval_call
(
 SEXP fn, SEXP args,
 double *t1, double *t2,
 double *x1, double *x2, int nvar,
 double *p, int npar,
 double *c, int ncov
 ) {

  SEXP var = args, ans, ob;
  int v;

  *(REAL(CAR(var))) = *t1; var = CDR(var);
  *(REAL(CAR(var))) = *t2; var = CDR(var);
  for (v = 0; v < nvar; v++, x1++, x2++) {
    *(REAL(CAR(var))) = *x1; var = CDR(var);
    *(REAL(CAR(var))) = *x2; var = CDR(var);
  }
  for (v = 0; v < npar; v++, p++, var=CDR(var)) *(REAL(CAR(var))) = *p;
  for (v = 0; v < ncov; v++, c++, var=CDR(var)) *(REAL(CAR(var))) = *c;

  PROTECT(ob = LCONS(fn,args));
  PROTECT(ans = eval(ob,CLOENV(fn)));

  UNPROTECT(2);
  return ans;

}

static R_INLINE SEXP ret_array (int nreps, int ntimes)
{
  int dim[2] = {nreps, ntimes};
  const char *dimnm[2] = {".id","time"};
  SEXP F;
  PROTECT(F = makearray(2,dim));
  fixdimnames(F,dimnm,2);
  UNPROTECT(1);
  return F;
}

// compute pdf of a sequence of elementary steps
static SEXP onestep_density
(
 SEXP func, SEXP x, SEXP times, SEXP params, SEXP covar,
 SEXP log, SEXP args, SEXP gnsi
 ) {

  pompfunmode mode = undef;
  int give_log;
  int nvars, npars, nrepsx, nrepsp, nreps, ntimes, ncovars;
  SEXP Snames, Pnames, Cnames;
  SEXP fn;
  SEXP F, cvec;
  double *cov;
  int *dim;

  ntimes = LENGTH(times);
  dim = INTEGER(GET_DIM(x)); nvars = dim[0]; nrepsx = dim[1];
  if (ntimes < 2)
    err("length(times) < 2: with no transitions, there is no work to do.");
  if (ntimes != dim[2])
    err("the length of 'times' and 3rd dimension of 'x' do not agree.");
  dim = INTEGER(GET_DIM(params)); npars = dim[0]; nrepsp = dim[1];

  give_log = *(INTEGER(AS_INTEGER(log)));

  // handle case with different numbers of states and parameters
  if (nrepsx != nrepsp && nrepsx % nrepsp != 0 && nrepsp % nrepsx != 0) {
    err("the larger number of replicates is not a multiple of smaller.");
  } else {
    nreps = (nrepsx > nrepsp) ? nrepsx : nrepsp;
  }

  PROTECT(Snames = GET_ROWNAMES(GET_DIMNAMES(x)));
  PROTECT(Pnames = GET_ROWNAMES(GET_DIMNAMES(params)));
  PROTECT(Cnames = get_covariate_names(covar));

  PROTECT(F = ret_array(nreps,ntimes-1));

  // set up the covariate table
  lookup_table_t covariate_table = make_covariate_table(covar,&ncovars);
  PROTECT(cvec = NEW_NUMERIC(ncovars));
  cov = REAL(cvec);

  PROTECT(fn = pomp_fun_handler(func,gnsi,&mode,Snames,Pnames,NA_STRING,Cnames));

  int nprotect = 6;
  double *t1 = REAL(times), *t2 = REAL(times)+1;
  double *x1p = REAL(x);
  double *x2p = REAL(x)+nrepsx*nvars;
  double *ft = REAL(F);

  switch (mode) {

  case Rfun: {

    PROTECT(args = add_args(args,Snames,Pnames,Cnames)); nprotect++;

    for (int k = 0; k < ntimes-1; k++, t1++, t2++) { // loop over times

      R_CheckUserInterrupt();

      // interpolate the covariates at time t1
      table_lookup(&covariate_table,*t1,cov);

      for (int j = 0; j < nreps; j++, ft++) {

        double *p = REAL(params)+npars*(j%nrepsp);
        double *x1 = x1p+nvars*(j%nrepsx);
        double *x2 = x2p+nvars*(j%nrepsx);

        *ft = *REAL(AS_NUMERIC(eval_call(fn,args,t1,t2,x1,x2,nvars,p,npars,cov,ncovars)));

        if (!give_log) *ft = exp(*ft);

      }

      x1p = x2p;
      x2p += nrepsx*nvars;

    }

  }

    break;

  case native: case regNative: {

    int *sidx, *pidx, *cidx;
    pomp_dprocess *ff = NULL;

    sidx = INTEGER(GET_SLOT(func,install("stateindex")));
    pidx = INTEGER(GET_SLOT(func,install("paramindex")));
    cidx = INTEGER(GET_SLOT(func,install("covarindex")));

    *((void **) (&ff)) = R_ExternalPtrAddr(fn);

    for (int k = 0; k < ntimes-1; k++, t1++, t2++) {

      R_CheckUserInterrupt();

      // interpolate the covariates at time t1
      table_lookup(&covariate_table,*t1,cov);

      for (int j = 0; j < nreps; j++, ft++) {

        double *p = REAL(params)+npars*(j%nrepsp);
        double *x1 = x1p+nvars*(j%nrepsx);
        double *x2 = x2p+nvars*(j%nrepsx);

        (*ff)(ft,x1,x2,*t1,*t2,p,sidx,pidx,cidx,cov);

        if (!give_log) *ft = exp(*ft);

      }

      x1p = x2p;
      x2p += nrepsx*nvars;

    }

  }

    break;

  default: {
    double *ft = REAL(F);
    int j, k;

    for (k = 0; k < ntimes-1; k++) { // loop over times
      for (j = 0; j < nreps; j++, ft++) { // loop over replicates
        *ft = R_NaReal;
      }
    }

    warn("'dprocess' unspecified: likelihood undefined.");

  }

  }

  UNPROTECT(nprotect);
  return F;
}

SEXP do_dprocess (SEXP object, SEXP x, SEXP times, SEXP params, SEXP log, SEXP gnsi)
{
  SEXP X, fn, args, covar;
  PROTECT(times=AS_NUMERIC(times));
  PROTECT(x = as_state_array(x));
  PROTECT(params = as_matrix(params));
  // extract the process function
  PROTECT(fn = GET_SLOT(object,install("dprocess")));
  // extract other arguments
  PROTECT(args = GET_SLOT(object,install("userdata")));
  PROTECT(covar = GET_SLOT(object,install("covar")));
  // evaluate the density
  PROTECT(X = onestep_density(fn,x,times,params,covar,log,args,gnsi));
  UNPROTECT(7);
  return X;
}
