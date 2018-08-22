// dear emacs, please treat this as -*- C++ -*-

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Arith.h>

#include "pomp_internal.h"

SEXP do_dprocess (SEXP object, SEXP x, SEXP times, SEXP params, SEXP log, SEXP gnsi)
{
  int nprotect = 0;
  int *xdim, npars, nvars, nreps, nrepsx, ntimes;
  SEXP X, fn, args, covar;

  PROTECT(gnsi = duplicate(gnsi)); nprotect++;

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

  {
    const char *dimnms[2] = {"rep","time"};
    fixdimnames(X,dimnms,2);
  }

  UNPROTECT(nprotect);
  return X;
}

// compute pdf of a sequence of elementary steps
SEXP onestep_density (SEXP func, SEXP x, SEXP times, SEXP params, SEXP covar,
  SEXP log, SEXP args, SEXP gnsi)
{
  int nprotect = 0;
  pompfunmode mode = undef;
  int give_log;
  int nvars, npars, nreps, ntimes, ncovars;
  pomp_onestep_pdf *ff = NULL;
  SEXP cvec, pvec = R_NilValue;
  SEXP t1vec = R_NilValue, t2vec = R_NilValue;
  SEXP x1vec = R_NilValue, x2vec = R_NilValue;
  SEXP Snames, Pnames, Cnames;
  SEXP fn, rho = R_NilValue, fcall = R_NilValue;
  SEXP F;
  int *pidx = 0, *sidx = 0, *cidx = 0;

  {
    int *dim;
    dim = INTEGER(GET_DIM(x)); nvars = dim[0]; nreps = dim[1];
    dim = INTEGER(GET_DIM(params)); npars = dim[0];
    ntimes = LENGTH(times);
  }

  PROTECT(Snames = GET_ROWNAMES(GET_DIMNAMES(x))); nprotect++;
  PROTECT(Pnames = GET_ROWNAMES(GET_DIMNAMES(params))); nprotect++;
  PROTECT(Cnames = get_covariate_names(covar)); nprotect++;

  // set up the covariate table
  lookup_table_t covariate_table = make_covariate_table(covar,&ncovars);
  // vector for interpolated covariates
  PROTECT(cvec = NEW_NUMERIC(ncovars)); nprotect++;
  SET_NAMES(cvec,Cnames);

  PROTECT(fn = pomp_fun_handler(func,gnsi,&mode)); nprotect++;

  give_log = *(INTEGER(log));

  switch (mode) {

  case Rfun:			// R function

    PROTECT(t1vec = NEW_NUMERIC(1)); nprotect++;
    PROTECT(t2vec = NEW_NUMERIC(1)); nprotect++;
    PROTECT(x1vec = NEW_NUMERIC(nvars)); nprotect++;
    SET_NAMES(x1vec,Snames);
    PROTECT(x2vec = NEW_NUMERIC(nvars)); nprotect++;
    SET_NAMES(x2vec,Snames);
    PROTECT(pvec = NEW_NUMERIC(npars)); nprotect++;
    SET_NAMES(pvec,Pnames);

    // set up the function call
    PROTECT(fcall = LCONS(cvec,args)); nprotect++;
    SET_TAG(fcall,install("covars"));
    PROTECT(fcall = LCONS(pvec,fcall)); nprotect++;
    SET_TAG(fcall,install("params"));
    PROTECT(fcall = LCONS(t2vec,fcall)); nprotect++;
    SET_TAG(fcall,install("t2"));
    PROTECT(fcall = LCONS(t1vec,fcall)); nprotect++;
    SET_TAG(fcall,install("t1"));
    PROTECT(fcall = LCONS(x2vec,fcall)); nprotect++;
    SET_TAG(fcall,install("x2"));
    PROTECT(fcall = LCONS(x1vec,fcall)); nprotect++;
    SET_TAG(fcall,install("x1"));
    PROTECT(fcall = LCONS(fn,fcall)); nprotect++;

    PROTECT(rho = (CLOENV(fn))); nprotect++;

    break;

  case native:			// native code

    // construct state, parameter, covariate indices
    sidx = INTEGER(PROTECT(matchnames(Snames,GET_SLOT(func,install("statenames")),"state variables"))); nprotect++;
    pidx = INTEGER(PROTECT(matchnames(Pnames,GET_SLOT(func,install("paramnames")),"parameters"))); nprotect++;
    cidx = INTEGER(PROTECT(matchnames(Cnames,GET_SLOT(func,install("covarnames")),"covariates"))); nprotect++;

    *((void **) (&ff)) = R_ExternalPtrAddr(fn);

    break;

  default:

    break;

  }

  // create array to hold results
  {
    int dim[2] = {nreps, ntimes-1};
    PROTECT(F = makearray(2,dim)); nprotect++;
  }

  switch (mode) {

  case Rfun:			// R function

  {
    double *cp = REAL(cvec);
    double *t1p = REAL(t1vec);
    double *t2p = REAL(t2vec);
    double *x1p = REAL(x1vec);
    double *x2p = REAL(x2vec);
    double *pp = REAL(pvec);
    double *t1s = REAL(times);
    double *t2s = t1s+1;
    double *x1s = REAL(x);
    double *x2s = x1s + nvars*nreps;
    double *ps;
    double *fs = REAL(F);
    int j, k;

    for (k = 0; k < ntimes-1; k++, t1s++, t2s++) { // loop over times

      R_CheckUserInterrupt();

      *t1p = *t1s; *t2p = *t2s;

      // interpolate the covariates at time t1, store the results in cvec
      table_lookup(&covariate_table,*t1p,cp);

      for (j = 0, ps = REAL(params); j < nreps; j++, fs++, x1s += nvars, x2s += nvars, ps += npars) { // loop over replicates

        memcpy(x1p,x1s,nvars*sizeof(double));
        memcpy(x2p,x2s,nvars*sizeof(double));
        memcpy(pp,ps,npars*sizeof(double));

        *fs = *(REAL(AS_NUMERIC(PROTECT(eval(fcall,rho)))));
        UNPROTECT(1);

        if (!give_log) *fs = exp(*fs);

      }
    }
  }

    break;

  case native:			// native code

    set_pomp_userdata(args);

    {
      double *t1s = REAL(times);
      double *t2s = t1s+1;
      double *x1s = REAL(x);
      double *x2s = x1s + nvars*nreps;
      double *fs = REAL(F);
      double *cp = REAL(cvec);
      double *ps;
      int j, k;

      for (k = 0; k < ntimes-1; k++, t1s++, t2s++) { // loop over times

        R_CheckUserInterrupt();

        // interpolate the covariates at time t1, store the results in cvec
        table_lookup(&covariate_table,*t1s,cp);

        for (j = 0, ps = REAL(params); j < nreps; j++, fs++, x1s += nvars, x2s += nvars, ps += npars) { // loop over replicates

          (*ff)(fs,x1s,x2s,*t1s,*t2s,ps,sidx,pidx,cidx,ncovars,cp);

          if (!give_log) *fs = exp(*fs);

        }
      }
    }

    unset_pomp_userdata();

    break;

  default:

    {
      double *fs = REAL(F);
      int j, k;

      for (k = 0; k < ntimes-1; k++) { // loop over times
        for (j = 0; j < nreps; j++, fs++) { // loop over replicates
          *fs = R_NaReal;
        }
      }
    }

  break;

  }

  UNPROTECT(nprotect);
  return F;
}
