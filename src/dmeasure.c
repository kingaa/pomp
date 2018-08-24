// dear emacs, please treat this as -*- C++ -*-

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Arith.h>

#include "pomp_internal.h"

SEXP do_dmeasure (SEXP object, SEXP y, SEXP x, SEXP times, SEXP params, SEXP log, SEXP gnsi)
{
  int nprotect = 0;
  pompfunmode mode = undef;
  int give_log;
  int ntimes, nvars, npars, ncovars, nreps, nrepsx, nrepsp, nobs;
  SEXP Snames, Pnames, Cnames, Onames;
  SEXP pompfun;
  SEXP cvec;
  SEXP fn, ans, fcall, rho = R_NilValue;
  SEXP F;
  int *sidx = 0, *pidx = 0, *cidx = 0, *oidx = 0;
  int *dim;
  lookup_table_t covariate_table;
  pomp_measure_model_density *ff = NULL;
  SEXP var;
  int v;

  PROTECT(times = AS_NUMERIC(times)); nprotect++;
  ntimes = length(times);
  if (ntimes < 1)
    errorcall(R_NilValue,"length('times') = 0, no work to do");

  PROTECT(y = as_matrix(y)); nprotect++;
  dim = INTEGER(GET_DIM(y));
  nobs = dim[0];

  if (ntimes != dim[1])
    errorcall(R_NilValue,"length of 'times' and 2nd dimension of 'y' do not agree");

  PROTECT(x = as_state_array(x)); nprotect++;
  dim = INTEGER(GET_DIM(x));
  nvars = dim[0]; nrepsx = dim[1];

  if (ntimes != dim[2])
    errorcall(R_NilValue,"length of 'times' and 3rd dimension of 'x' do not agree");

  PROTECT(params = as_matrix(params)); nprotect++;
  dim = INTEGER(GET_DIM(params));
  npars = dim[0]; nrepsp = dim[1];

  nreps = (nrepsp > nrepsx) ? nrepsp : nrepsx;

  if ((nreps % nrepsp != 0) || (nreps % nrepsx != 0))
    errorcall(R_NilValue,"larger number of replicates is not a multiple of smaller");

  PROTECT(Onames = GET_ROWNAMES(GET_DIMNAMES(y))); nprotect++;
  PROTECT(Snames = GET_ROWNAMES(GET_DIMNAMES(x))); nprotect++;
  PROTECT(Pnames = GET_ROWNAMES(GET_DIMNAMES(params))); nprotect++;
  PROTECT(Cnames = get_covariate_names(GET_SLOT(object,install("covar")))); nprotect++;

  give_log = *(INTEGER(AS_INTEGER(log)));

  // set up the covariate table
  covariate_table = make_covariate_table(GET_SLOT(object,install("covar")),&ncovars);

  // vector for interpolated covariates
  PROTECT(cvec = NEW_NUMERIC(ncovars)); nprotect++;
  SET_NAMES(cvec,Cnames);

  // extract the user-defined function
  PROTECT(pompfun = GET_SLOT(object,install("dmeasure"))); nprotect++;
  PROTECT(fn = pomp_fun_handler(pompfun,gnsi,&mode)); nprotect++;

  // extract 'userdata' as pairlist
  PROTECT(fcall = VectorToPairList(GET_SLOT(object,install("userdata")))); nprotect++;

  // first do setup
  switch (mode) {

  case Rfun:			// R function

    PROTECT(fcall = LCONS(log,fcall)); nprotect++;
    SET_TAG(fcall,install("log"));
    for (v = LENGTH(Cnames)-1; v >= 0; v--) {
      PROTECT(var = NEW_NUMERIC(1)); nprotect++;
      PROTECT(fcall = LCONS(var,fcall)); nprotect++;
      SET_TAG(fcall,install(CHAR(STRING_ELT(Cnames,v))));
    }
    for (v = LENGTH(Pnames)-1; v >= 0; v--) {
      PROTECT(var = NEW_NUMERIC(1)); nprotect++;
      PROTECT(fcall = LCONS(var,fcall)); nprotect++;
      SET_TAG(fcall,install(CHAR(STRING_ELT(Pnames,v))));
    }
    for (v = LENGTH(Snames)-1; v >= 0; v--) {
      PROTECT(var = NEW_NUMERIC(1)); nprotect++;
      PROTECT(fcall = LCONS(var,fcall)); nprotect++;
      SET_TAG(fcall,install(CHAR(STRING_ELT(Snames,v))));
    }
    for (v = LENGTH(Onames)-1; v >= 0; v--) {
      PROTECT(var = NEW_NUMERIC(1)); nprotect++;
      PROTECT(fcall = LCONS(var,fcall)); nprotect++;
      SET_TAG(fcall,install(CHAR(STRING_ELT(Onames,v))));
    }
    PROTECT(var = NEW_NUMERIC(1)); nprotect++;
    PROTECT(fcall = LCONS(var,fcall)); nprotect++;
    SET_TAG(fcall,install("t"));

    PROTECT(fcall = LCONS(fn,fcall)); nprotect++;

    // get the function's environment
    PROTECT(rho = (CLOENV(fn))); nprotect++;

    break;

  case native:			// native code

    // construct state, parameter, covariate, observable indices
    oidx = INTEGER(PROTECT(name_index(Onames,pompfun,"obsnames","observables"))); nprotect++;
    sidx = INTEGER(PROTECT(name_index(Snames,pompfun,"statenames","state variables"))); nprotect++;
    pidx = INTEGER(PROTECT(name_index(Pnames,pompfun,"paramnames","parameters"))); nprotect++;
    cidx = INTEGER(PROTECT(name_index(Cnames,pompfun,"covarnames","covariates"))); nprotect++;

    // address of native routine
    *((void **) (&ff)) = R_ExternalPtrAddr(fn);

    break;

  default:

    break;

  }

  // create array to store results
  {
    int dim[2] = {nreps, ntimes};
    const char *dimnm[2] = {"rep","time"};
    PROTECT(F = makearray(2,dim)); nprotect++;
    fixdimnames(F,dimnm,2);
  }

  // now do computations
  switch (mode) {

  case Rfun:			// R function

  {
    int first = 1;
    double *ys = REAL(y);
    double *xs = REAL(x);
    double *ps = REAL(params);
    double *cp = REAL(cvec);
    double *ft = REAL(F);
    double *time = REAL(times);
    int j, k;

    for (k = 0; k < ntimes; k++, time++, ys += nobs) { // loop over times

      R_CheckUserInterrupt();	// check for user interrupt

      table_lookup(&covariate_table,*time,cp); // interpolate the covariates

      for (j = 0; j < nreps; j++, ft++) { // loop over replicates

        var = CDR(fcall);
        // copy the states and parameters into place
        *(REAL(CAR(var))) = *time; var = CDR(var);
        for (v = 0; v < nobs; v++, var=CDR(var))
          *(REAL(CAR(var))) = ys[v];
        for (v = 0; v < nvars; v++, var=CDR(var))
          *(REAL(CAR(var))) = xs[v+nvars*((j%nrepsx)+nrepsx*k)];
        for (v = 0; v < npars; v++, var=CDR(var))
          *(REAL(CAR(var))) = ps[v+npars*(j%nrepsp)];
        for (v = 0; v < ncovars; v++, var=CDR(var))
          *(REAL(CAR(var))) = cp[v];

        if (first) {
          // evaluate the call
          PROTECT(ans = eval(fcall,rho)); nprotect++;
          if (LENGTH(ans) != 1)
            errorcall(R_NilValue,"user 'dmeasure' returns a vector of length %d when it should return a scalar",LENGTH(ans));

          *ft = *(REAL(AS_NUMERIC(ans)));

          first = 0;

        } else {

          *ft = *(REAL(AS_NUMERIC(PROTECT(eval(fcall,rho)))));
          UNPROTECT(1);
        }

      }
    }
  }

    break;

  case native:			// native code

    set_pomp_userdata(fcall);

    {
      double *yp = REAL(y);
      double *xs = REAL(x);
      double *ps = REAL(params);
      double *cp = REAL(cvec);
      double *ft = REAL(F);
      double *time = REAL(times);
      double *xp, *pp;
      int j, k;

      for (k = 0; k < ntimes; k++, time++, yp += nobs) { // loop over times

        R_CheckUserInterrupt();	// check for user interrupt

        // interpolate the covar functions for the covariates
        table_lookup(&covariate_table,*time,cp);

        for (j = 0; j < nreps; j++, ft++) { // loop over replicates

          xp = &xs[nvars*((j%nrepsx)+nrepsx*k)];
          pp = &ps[npars*(j%nrepsp)];

          (*ff)(ft,yp,xp,pp,give_log,oidx,sidx,pidx,cidx,ncovars,cp,*time);

        }
      }
    }

    unset_pomp_userdata();

    break;

  default:

    {
      double *ft = REAL(F);
      int j, k;

      for (k = 0; k < ntimes; k++) { // loop over times

        for (j = 0; j < nreps; j++, ft++) { // loop over replicates

          *ft = R_NaReal;

        }
      }
    }

  break;

  }

  UNPROTECT(nprotect);
  return F;
}
