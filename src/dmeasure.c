// dear emacs, please treat this as -*- C++ -*-

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Arith.h>

#include "pomp_internal.h"

static R_INLINE SEXP dmeas_args (SEXP args, SEXP Onames, SEXP Snames,
  SEXP Pnames, SEXP Cnames, SEXP log)
{
  int nprotect = 0;
  SEXP var;
  int v;

  // we construct the call from end to beginning
  // 'log', covariates, parameter, states, observables, then time

  // 'log' is a needed argument
  PROTECT(args = LCONS(AS_LOGICAL(log),args)); nprotect++;
  SET_TAG(args,install("log"));

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

static R_INLINE SEXP eval_call (
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

static R_INLINE SEXP dmeas_array (int nreps, int ntimes) {
  int dim[2] = {nreps, ntimes};
  const char *dimnm[2] = {"rep","time"};
  SEXP F;
  PROTECT(F = makearray(2,dim));
  fixdimnames(F,dimnm,2);
  UNPROTECT(1);
  return F;
}

SEXP do_dmeasure (SEXP object, SEXP y, SEXP x, SEXP times, SEXP params, SEXP log, SEXP gnsi)
{
  int nprotect = 0;
  pompfunmode mode = undef;
  int ntimes, nvars, npars, ncovars, nreps, nrepsx, nrepsp, nobs;
  SEXP Snames, Pnames, Cnames, Onames;
  SEXP cvec, pompfun;
  SEXP fn, args, ans;
  SEXP F;
  int *dim;
  lookup_table_t covariate_table;
  double *cov;

  PROTECT(times = AS_NUMERIC(times)); nprotect++;
  ntimes = length(times);
  if (ntimes < 1) errorcall(R_NilValue,"length('times') = 0, no work to do.");

  PROTECT(y = as_matrix(y)); nprotect++;
  dim = INTEGER(GET_DIM(y));
  nobs = dim[0];

  if (ntimes != dim[1]) errorcall(R_NilValue,"length of 'times' and 2nd dimension of 'y' do not agree.");

  PROTECT(x = as_state_array(x)); nprotect++;
  dim = INTEGER(GET_DIM(x));
  nvars = dim[0]; nrepsx = dim[1];

  if (ntimes != dim[2])
    errorcall(R_NilValue,"length of 'times' and 3rd dimension of 'x' do not agree.");

  PROTECT(params = as_matrix(params)); nprotect++;
  dim = INTEGER(GET_DIM(params));
  npars = dim[0]; nrepsp = dim[1];

  nreps = (nrepsp > nrepsx) ? nrepsp : nrepsx;

  if ((nreps % nrepsp != 0) || (nreps % nrepsx != 0))
    errorcall(R_NilValue,"larger number of replicates is not a multiple of smaller.");

  PROTECT(Onames = GET_ROWNAMES(GET_DIMNAMES(y))); nprotect++;
  PROTECT(Snames = GET_ROWNAMES(GET_DIMNAMES(x))); nprotect++;
  PROTECT(Pnames = GET_ROWNAMES(GET_DIMNAMES(params))); nprotect++;
  PROTECT(Cnames = get_covariate_names(GET_SLOT(object,install("covar")))); nprotect++;

  // set up the covariate table
  covariate_table = make_covariate_table(GET_SLOT(object,install("covar")),&ncovars);
  PROTECT(cvec = NEW_NUMERIC(ncovars)); nprotect++;
  cov = REAL(cvec);

  // extract the user-defined function
  PROTECT(pompfun = GET_SLOT(object,install("dmeasure"))); nprotect++;
  PROTECT(fn = pomp_fun_handler(pompfun,gnsi,&mode,Snames,Pnames,Onames,Cnames)); nprotect++;

  // extract 'userdata' as pairlist
  PROTECT(args = VectorToPairList(GET_SLOT(object,install("userdata")))); nprotect++;

  // create array to store results
  PROTECT(F = dmeas_array(nreps,ntimes)); nprotect++;

  switch (mode) {

  case Rfun: {
    double *ys = REAL(y), *xs = REAL(x), *ps = REAL(params), *time = REAL(times);
    double *ft = REAL(F);
    int j, k;

    // build argument list
    PROTECT(args = dmeas_args(args,Onames,Snames,Pnames,Cnames,log)); nprotect++;

    for (k = 0; k < ntimes; k++, time++, ys += nobs) { // loop over times

      R_CheckUserInterrupt();	// check for user interrupt

      table_lookup(&covariate_table,*time,cov); // interpolate the covariates

      for (j = 0; j < nreps; j++, ft++) { // loop over replicates

        // evaluate the call
        PROTECT(
          ans = eval_call(
            fn,args,
            time,
            ys,nobs,
            xs+nvars*((j%nrepsx)+nrepsx*k),nvars,
            ps+npars*(j%nrepsp),npars,
            cov,ncovars
          )
        );

        if (k == 0 && j == 0 && LENGTH(ans) != 1)
          errorcall(R_NilValue,"user 'dmeasure' returns a vector of length %d when it should return a scalar.",LENGTH(ans));

        *ft = *(REAL(AS_NUMERIC(ans)));

        UNPROTECT(1);

      }
    }
  }

    break;

  case native: case regNative: {
    int *oidx, *sidx, *pidx, *cidx;
    int give_log;
    pomp_measure_model_density *ff = NULL;
    double *yp = REAL(y), *xs = REAL(x), *ps = REAL(params), *time = REAL(times);
    double *ft = REAL(F);
    double *xp, *pp;
    int j, k;

    // extract state, parameter, covariate, observable indices
    sidx = INTEGER(GET_SLOT(pompfun,install("stateindex")));
    pidx = INTEGER(GET_SLOT(pompfun,install("paramindex")));
    oidx = INTEGER(GET_SLOT(pompfun,install("obsindex")));
    cidx = INTEGER(GET_SLOT(pompfun,install("covarindex")));

    give_log = *(INTEGER(AS_INTEGER(log)));

    // address of native routine
    *((void **) (&ff)) = R_ExternalPtrAddr(fn);

    set_pomp_userdata(args);

    for (k = 0; k < ntimes; k++, time++, yp += nobs) { // loop over times

      R_CheckUserInterrupt();	// check for user interrupt

      // interpolate the covar functions for the covariates
      table_lookup(&covariate_table,*time,cov);

      for (j = 0; j < nreps; j++, ft++) { // loop over replicates

        xp = &xs[nvars*((j%nrepsx)+nrepsx*k)];
        pp = &ps[npars*(j%nrepsp)];

        (*ff)(ft,yp,xp,pp,give_log,oidx,sidx,pidx,cidx,cov,*time);

      }
    }

    unset_pomp_userdata();

  }

    break;

  default: {
    double *ft = REAL(F);
    int j, k;

    for (k = 0; k < ntimes; k++) { // loop over times

      for (j = 0; j < nreps; j++, ft++) { // loop over replicates

        *ft = R_NaReal;

      }
    }
  }

  }

  UNPROTECT(nprotect);
  return F;
}
