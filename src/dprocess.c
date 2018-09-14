// dear emacs, please treat this as -*- C++ -*-

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Arith.h>

#include "pomp_internal.h"

static R_INLINE SEXP add_args (SEXP args, SEXP Snames, SEXP Pnames, SEXP Cnames)
{
  int nprotect = 0;
  SEXP S1names, S2names;
  SEXP var;
  int v;

  PROTECT(S1names = paste0(Snames,mkString("_1"))); nprotect++;
  PROTECT(S2names = paste0(Snames,mkString("_2"))); nprotect++;

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
    SET_TAG(args,install(CHAR(STRING_ELT(S2names,v))));
    PROTECT(var = NEW_NUMERIC(1)); nprotect++;
    PROTECT(args = LCONS(var,args)); nprotect++;
    SET_TAG(args,install(CHAR(STRING_ELT(S1names,v))));
  }

  // Time
  PROTECT(var = NEW_NUMERIC(1)); nprotect++;
  PROTECT(args = LCONS(var,args)); nprotect++;
  SET_TAG(args,install("t_2"));
  PROTECT(var = NEW_NUMERIC(1)); nprotect++;
  PROTECT(args = LCONS(var,args)); nprotect++;
  SET_TAG(args,install("t_1"));

  UNPROTECT(nprotect);
  return args;

}

static R_INLINE SEXP eval_call (
    SEXP fn, SEXP args,
    double *t1, double *t2,
    double *x1, double *x2, int nvar,
    double *p, int npar,
    double *c, int ncov)
{

  SEXP var = args, ans;
  int v;

  *(REAL(CAR(var))) = *t1; var = CDR(var);
  *(REAL(CAR(var))) = *t2; var = CDR(var);
  for (v = 0; v < nvar; v++, x1++, x2++) {
    *(REAL(CAR(var))) = *x1; var = CDR(var);
    *(REAL(CAR(var))) = *x2; var = CDR(var);
  }
  for (v = 0; v < npar; v++, p++, var=CDR(var)) *(REAL(CAR(var))) = *p;
  for (v = 0; v < ncov; v++, c++, var=CDR(var)) *(REAL(CAR(var))) = *c;

  PROTECT(ans = eval(LCONS(fn,args),CLOENV(fn)));

  UNPROTECT(1);
  return ans;

}

static R_INLINE SEXP ret_array (int nreps, int ntimes)
{
  int dim[2] = {nreps, ntimes};
  const char *dimnm[2] = {"rep","time"};
  SEXP F;
  PROTECT(F = makearray(2,dim));
  fixdimnames(F,dimnm,2);
  UNPROTECT(1);
  return F;
}

// compute pdf of a sequence of elementary steps
static SEXP onestep_density (SEXP func, SEXP x, SEXP times, SEXP params, SEXP covar,
  SEXP log, SEXP args, SEXP gnsi)
{
  int nprotect = 0;
  pompfunmode mode = undef;
  int give_log;
  int nvars, npars, nreps, ntimes, ncovars;
  SEXP Snames, Pnames, Cnames;
  SEXP fn;
  SEXP F, cvec;
  double *cov;
  int *dim;

  dim = INTEGER(GET_DIM(x)); nvars = dim[0]; nreps = dim[1];
  dim = INTEGER(GET_DIM(params)); npars = dim[0];
  ntimes = LENGTH(times);

  PROTECT(Snames = GET_ROWNAMES(GET_DIMNAMES(x))); nprotect++;
  PROTECT(Pnames = GET_ROWNAMES(GET_DIMNAMES(params))); nprotect++;
  PROTECT(Cnames = get_covariate_names(covar)); nprotect++;

  PROTECT(F = ret_array(nreps,ntimes-1)); nprotect++;

  // set up the covariate table
  lookup_table_t covariate_table = make_covariate_table(covar,&ncovars);
  PROTECT(cvec = NEW_NUMERIC(ncovars)); nprotect++;
  cov = REAL(cvec);

  PROTECT(fn = pomp_fun_handler(func,gnsi,&mode,Snames,Pnames,NA_STRING,Cnames)); nprotect++;

  give_log = *(INTEGER(log));

  switch (mode) {

  case Rfun: {

    SEXP ans;
    double *t1 = REAL(times), *t2 = REAL(times)+1;
    double *ps;
    double *x1 = REAL(x), *x2 = REAL(x)+nvars*nreps;
    double *ft = REAL(F);
    int j, k;

    PROTECT(args = add_args(args,Snames,Pnames,Cnames)); nprotect++;

    for (k = 0; k < ntimes-1; k++, t1++, t2++) { // loop over times

      R_CheckUserInterrupt();

      // interpolate the covariates at time t1
      table_lookup(&covariate_table,*t1,cov);

      for (j = 0, ps = REAL(params);
        j < nreps;
        j++, ft++, x1 += nvars, x2 += nvars, ps += npars) {

        PROTECT(ans = eval_call(fn,args,t1,t2,x1,x2,nvars,ps,npars,cov,ncovars));

        *ft = *(REAL(AS_NUMERIC(ans)));

        UNPROTECT(1);

        if (!give_log) *ft = exp(*ft);

      }
    }


  }

    break;

  case native: case regNative: {

    int *sidx, *pidx, *cidx;
    double *t1 = REAL(times), *t2 = REAL(times)+1;
    double *ps = REAL(params);
    double *x1 = REAL(x), *x2 = REAL(x)+nvars*nreps;
    double *ft = REAL(F);
    pomp_onestep_pdf *ff = NULL;
    int j, k;

    sidx = INTEGER(GET_SLOT(func,install("stateindex")));
    pidx = INTEGER(GET_SLOT(func,install("paramindex")));
    cidx = INTEGER(GET_SLOT(func,install("covarindex")));

    *((void **) (&ff)) = R_ExternalPtrAddr(fn);

    set_pomp_userdata(args);

    for (k = 0; k < ntimes-1; k++, t1++, t2++) {

      R_CheckUserInterrupt();

      // interpolate the covariates at time t1
      table_lookup(&covariate_table,*t1,cov);

      for (j = 0, ps = REAL(params); j < nreps; j++, ft++, x1 += nvars, x2 += nvars, ps += npars) {

        (*ff)(ft,x1,x2,*t1,*t2,ps,sidx,pidx,cidx,cov);

        if (!give_log) *ft = exp(*ft);

      }
    }

    unset_pomp_userdata();
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

    warningcall(R_NilValue,"'dprocess' unspecified: likelihood undefined.");

  }

  }

  UNPROTECT(nprotect);
  return F;
}

SEXP do_dprocess (SEXP object, SEXP x, SEXP times, SEXP params, SEXP log, SEXP gnsi)
{
  int nprotect = 0;
  int *xdim, npars, nvars, nreps, nrepsx, ntimes;
  SEXP X, fn, args, covar;

  PROTECT(times=AS_NUMERIC(times)); nprotect++;
  ntimes = length(times);
  if (ntimes < 2)
    errorcall(R_NilValue,"length(times)<2: with no transitions, there is no work to do.");

  PROTECT(x = as_state_array(x)); nprotect++;
  xdim = INTEGER(GET_DIM(x));
  nvars = xdim[0]; nrepsx = xdim[1];
  if (ntimes != xdim[2])
    errorcall(R_NilValue,"the length of 'times' and 3rd dimension of 'x' do not agree.");

  PROTECT(params = as_matrix(params)); nprotect++;
  xdim = INTEGER(GET_DIM(params));
  npars = xdim[0]; nreps = xdim[1];

  if (nrepsx > nreps) {         // more states than parameters
    if (nrepsx % nreps != 0) {
      errorcall(R_NilValue,"the larger number of replicates is not a multiple of smaller.");
    } else {
      SEXP copy;
      double *src, *tgt;
      int dims[2];
      int j, k;
      dims[0] = npars; dims[1] = nrepsx;
      PROTECT(copy = duplicate(params)); nprotect++;
      PROTECT(params = makearray(2,dims)); nprotect++;
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
      errorcall(R_NilValue,"the larger number of replicates is not a multiple of smaller.");
    } else {
      SEXP copy;
      double *src, *tgt;
      int dims[3];
      int i, j, k;
      dims[0] = nvars; dims[1] = nreps; dims[2] = ntimes;
      PROTECT(copy = duplicate(x)); nprotect++;
      PROTECT(x = makearray(3,dims)); nprotect++;
      setrownames(x,GET_ROWNAMES(GET_DIMNAMES(copy)),3);
      src = REAL(copy);
      tgt = REAL(x);
      for (i = 0; i < ntimes; i++) {
        for (j = 0; j < nreps; j++) {
          for (k = 0; k < nvars; k++, tgt++) {
            *tgt = src[k+nvars*((j%nrepsx)+nrepsx*i)];
          }
        }
      }
    }
  }

  // extract the process function
  PROTECT(fn = GET_SLOT(object,install("dprocess"))); nprotect++;
  // extract other arguments
  PROTECT(args = VectorToPairList(GET_SLOT(object,install("userdata")))); nprotect++;
  PROTECT(covar = GET_SLOT(object,install("covar"))); nprotect++;

  PROTECT(X = onestep_density(fn,x,times,params,covar,log,args,gnsi)); nprotect++;

  UNPROTECT(nprotect);
  return X;
}
