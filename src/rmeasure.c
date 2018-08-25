// dear emacs, please treat this as -*- C++ -*-

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Arith.h>

#include "pomp_internal.h"

SEXP do_rmeasure (SEXP object, SEXP x, SEXP times, SEXP params, SEXP gnsi)
{
  int nprotect = 0;
  pompfunmode mode = undef;
  int ntimes, nvars, npars, ncovars, nreps, nrepsx, nrepsp;
  int nobs = 0;
  SEXP Snames, Pnames, Cnames, Onames = R_NilValue;
  SEXP cvec, indices;
  SEXP fn, args;
  SEXP pompfun;
  SEXP Y = R_NilValue;
  int *dim;
  int *sidx = 0, *pidx = 0, *cidx = 0, *oidx = 0;
  lookup_table_t covariate_table;
  pomp_measure_model_simulator *ff = NULL;

  PROTECT(times = AS_NUMERIC(times)); nprotect++;
  ntimes = length(times);
  if (ntimes < 1)
    errorcall(R_NilValue,"length('times') = 0, no work to do.");

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

  PROTECT(Snames = GET_ROWNAMES(GET_DIMNAMES(x))); nprotect++;
  PROTECT(Pnames = GET_ROWNAMES(GET_DIMNAMES(params))); nprotect++;
  PROTECT(Cnames = get_covariate_names(GET_SLOT(object,install("covar")))); nprotect++;

  // set up the covariate table
  covariate_table = make_covariate_table(GET_SLOT(object,install("covar")),&ncovars);

  // vector for interpolated covariates
  PROTECT(cvec = NEW_NUMERIC(ncovars)); nprotect++;
  SET_NAMES(cvec,Cnames);

  // extract the user-defined function
  PROTECT(pompfun = GET_SLOT(object,install("rmeasure"))); nprotect++;
  PROTECT(fn = pomp_fun_handler(pompfun,gnsi,&mode)); nprotect++;

  // extract 'userdata' as pairlist
  PROTECT(args = VectorToPairList(GET_SLOT(object,install("userdata")))); nprotect++;

  // first do setup
  switch (mode) {
  case Rfun:			// use R function

    PROTECT(args = pomp_fun_args(args,R_NilValue,Snames,Pnames,Cnames)); nprotect++;

    break;

  case native: case regNative:

    // construct observable, state, parameter covariate indices
    PROTECT(Onames = GET_SLOT(pompfun,install("obsnames"))); nprotect++;
    PROTECT(indices = pomp_fun_indices(pompfun,Onames,Snames,Pnames,Cnames)); nprotect++;
    oidx = INTEGER(VECTOR_ELT(indices,0));
    sidx = INTEGER(VECTOR_ELT(indices,1));
    pidx = INTEGER(VECTOR_ELT(indices,2));
    cidx = INTEGER(VECTOR_ELT(indices,3));
    nobs = LENGTH(Onames);

    // address of native routine
    *((void **) (&ff)) = R_ExternalPtrAddr(fn);

    break;

  default:

    PROTECT(Onames = GET_SLOT(pompfun,install("obsnames"))); nprotect++;
    nobs = LENGTH(Onames);

    break;

  }

  // now do computations
  switch (mode) {

  case Rfun:			// R function

  {
    double *yt = 0;
    double *time = REAL(times);
    double *cs = REAL(cvec);
    double *xs = REAL(x);
    double *ps = REAL(params);
    double *ys;
    SEXP ans;
    int i, j, k;

    for (k = 0; k < ntimes; k++, time++) { // loop over times

      R_CheckUserInterrupt();	// check for user interrupt

      table_lookup(&covariate_table,*time,cs); // interpolate the covariates

      for (j = 0; j < nreps; j++) { // loop over replicates

        if (k == 0 && j == 0) {

          const char *dimnm[3] = {"variable","rep","time"};
          int dim[3];

          PROTECT(
            ans = eval_pomp_fun_R_call(
              fn,args,
              time,
              0,0,
              xs+nvars*((j%nrepsx)+nrepsx*k),nvars,
              ps+npars*(j%nrepsp),npars,
              cs,ncovars
            )
          ); nprotect++;

          nobs = LENGTH(ans);

          PROTECT(Onames = GET_NAMES(ans)); nprotect++;
          if (isNull(Onames))
            errorcall(R_NilValue,"'rmeasure' must return a named numeric vector.");

          dim[0] = nobs; dim[1] = nreps; dim[2] = ntimes;
          PROTECT(Y = makearray(3,dim)); nprotect++;
          setrownames(Y,Onames,3);
          fixdimnames(Y,dimnm,3);

          yt = REAL(Y);
          ys = REAL(AS_NUMERIC(ans));

          memcpy(yt,ys,nobs*sizeof(double));
          yt += nobs;

        } else {

          PROTECT(
            ans = eval_pomp_fun_R_call(
              fn,args,
              time,
              0,0,
              xs+nvars*((j%nrepsx)+nrepsx*k),nvars,
              ps+npars*(j%nrepsp),npars,
              cs,ncovars
            )
          );

          if (LENGTH(ans) != nobs)
            errorcall(R_NilValue,"'rmeasure' returns variable-length results.");

          ys = REAL(AS_NUMERIC(ans));

          memcpy(yt,ys,nobs*sizeof(double));
          yt += nobs;

          UNPROTECT(1);

        }

      }
    }
  }

    break;

  case native: case regNative:			// native routine

  {
    int dim[3] = {nobs, nreps, ntimes};
    const char *dimnm[3] = {"variable","rep","time"};

    PROTECT(Y = makearray(3,dim)); nprotect++;
    setrownames(Y,Onames,3);
    fixdimnames(Y,dimnm,3);

    double *yt = REAL(Y);
    double *time = REAL(times);
    double *xs = REAL(x);
    double *ps = REAL(params);
    double *cp = REAL(cvec);
    double *xp, *pp;
    int j, k;

    set_pomp_userdata(args);
    GetRNGstate();

    for (k = 0; k < ntimes; k++, time++) { // loop over times

      R_CheckUserInterrupt();	// check for user interrupt

      // interpolate the covar functions for the covariates
      table_lookup(&covariate_table,*time,cp);

      for (j = 0; j < nreps; j++, yt += nobs) { // loop over replicates

        xp = &xs[nvars*((j%nrepsx)+nrepsx*k)];
        pp = &ps[npars*(j%nrepsp)];

        (*ff)(yt,xp,pp,oidx,sidx,pidx,cidx,ncovars,cp,*time);

      }
    }

    PutRNGstate();
    unset_pomp_userdata();
  }

    break;

  default:

  {
    int dim[3] = {nobs, nreps, ntimes};
    const char *dimnm[3] = {"variable","rep","time"};
    double *yt = 0;
    int i, n = nobs*nreps*ntimes;

    PROTECT(Y = makearray(3,dim)); nprotect++;
    setrownames(Y,Onames,3);
    fixdimnames(Y,dimnm,3);

    for (i = 0, yt = REAL(Y); i < n; i++, yt++) *yt = R_NaReal;

  }

  break;

  }

  UNPROTECT(nprotect);
  return Y;
}
