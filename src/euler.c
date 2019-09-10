// dear emacs, please treat this as -*- C++ -*-

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Constants.h>
#include <string.h>

#include "pomp_internal.h"

static R_INLINE SEXP add_args (SEXP args, SEXP Snames, SEXP Pnames, SEXP Cnames)
{

  SEXP var;
  int v;

  PROTECT(args);

  // we construct the call from end to beginning
  // delta.t, covariates, parameter, states, then time

  // 'delta.t'
  var = NEW_NUMERIC(1);
  args = LCONS(var,args);
  UNPROTECT(1);
  PROTECT(args);
  SET_TAG(args,install("delta.t"));

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
  SET_TAG(args,install("t"));

  UNPROTECT(1);
  return args;

}

static R_INLINE SEXP eval_call (
    SEXP fn, SEXP args,
    double *t, double *dt,
    double *x, int nvar,
    double *p, int npar,
    double *c, int ncov)
{

  SEXP var = args, ans, ob;
  int v;

  *(REAL(CAR(var))) = *t; var = CDR(var);
  for (v = 0; v < nvar; v++, x++, var=CDR(var)) *(REAL(CAR(var))) = *x;
  for (v = 0; v < npar; v++, p++, var=CDR(var)) *(REAL(CAR(var))) = *p;
  for (v = 0; v < ncov; v++, c++, var=CDR(var)) *(REAL(CAR(var))) = *c;
  *(REAL(CAR(var))) = *dt; var = CDR(var);

  PROTECT(ob = LCONS(fn,args));
  PROTECT(ans = eval(ob,CLOENV(fn)));

  UNPROTECT(2);
  return ans;

}

static R_INLINE SEXP ret_array (int n, int nreps, int ntimes, SEXP names)
{
  int dim[3] = {n, nreps, ntimes};
  const char *dimnm[3] = {"variable", "rep", "time"};
  SEXP Y;

  PROTECT(Y = makearray(3,dim));
  setrownames(Y,names,3);
  fixdimnames(Y,dimnm,3);

  UNPROTECT(1);
  return Y;

}

SEXP euler_model_simulator (SEXP func, SEXP xstart, SEXP tstart, SEXP times, SEXP params,
  double deltat, rprocmode method, SEXP accumvars, SEXP covar, SEXP args, SEXP gnsi)
{

  pompfunmode mode = undef;
  int nvars, npars, nreps, ntimes, nzeros, ncovars;
  double *cov, t0;
  SEXP cvec, X, fn;

  if (deltat <= 0) errorcall(R_NilValue,"'delta.t' should be a positive number."); // #nocov

  int *dim;
  dim = INTEGER(GET_DIM(xstart)); nvars = dim[0]; nreps = dim[1];
  dim = INTEGER(GET_DIM(params)); npars = dim[0];
  ntimes = LENGTH(times);

  PROTECT(tstart = AS_NUMERIC(tstart));
  PROTECT(times = AS_NUMERIC(times));
  t0 = *(REAL(tstart));
  if (t0 > *(REAL(times))) errorcall(R_NilValue,"'t0' must be no later than 'times[1]'.");

  SEXP Snames, Pnames, Cnames;
  PROTECT(Snames = GET_ROWNAMES(GET_DIMNAMES(xstart)));
  PROTECT(Pnames = GET_ROWNAMES(GET_DIMNAMES(params)));
  PROTECT(Cnames = get_covariate_names(covar));

  // set up the covariate table
  lookup_table_t covariate_table = make_covariate_table(covar,&ncovars);
  PROTECT(cvec = NEW_NUMERIC(ncovars));
  cov = REAL(cvec);

  // indices of accumulator variables
  nzeros = LENGTH(accumvars);
  int *zidx = INTEGER(PROTECT(matchnames(Snames,accumvars,"state variables")));

  // extract user function
  PROTECT(fn = pomp_fun_handler(func,gnsi,&mode,Snames,Pnames,NA_STRING,Cnames));

  // array to hold results
  PROTECT(X = ret_array(nvars,nreps,ntimes,Snames));

  // copy the start values into the result array
  memcpy(REAL(X),REAL(xstart),nvars*nreps*sizeof(double));

  // set up

  int nprotect = 9;
  int *pidx = 0, *sidx = 0, *cidx = 0;
  pomp_onestep_sim *ff = NULL;

  switch (mode) {

  case Rfun: {

    // construct list of all arguments
    PROTECT(args = add_args(args,Snames,Pnames,Cnames)); nprotect++;

  }

    break;

  case native: case regNative: {

    // construct state, parameter, covariate indices
    sidx = INTEGER(GET_SLOT(func,install("stateindex")));
    pidx = INTEGER(GET_SLOT(func,install("paramindex")));
    cidx = INTEGER(GET_SLOT(func,install("covarindex")));

    *((void **) (&ff)) = R_ExternalPtrAddr(fn);

    set_pomp_userdata(args);
    GetRNGstate();

  }

    break;

  default: // # nocov

    errorcall(R_NilValue,"unrecognized 'mode' %d",mode); // # nocov

  }

  // main computation loop
  int step;
  double *xt, *time, t;
  int first = 1;

  for (step = 0, xt = REAL(X), time = REAL(times), t = t0;
    step < ntimes;
    step++, xt += nvars*nreps) {

    double dt;
    int nstep = 0;
    int i, j, k;

    R_CheckUserInterrupt();

    if (t > time[step]) errorcall(R_NilValue,"'times' must be an increasing sequence.");  // #nocov

    // set accumulator variables to zero
    for (j = 0; j < nreps; j++)
      for (i = 0; i < nzeros; i++)
        xt[zidx[i]+nvars*j] = 0.0;

    // determine size and number of time-steps
    switch (method) {
    case onestep: default:	// one step
      dt = time[step]-t;
      nstep = (dt > 0) ? 1 : 0;
      break;
    case discrete:			// fixed step
      dt = deltat;
      nstep = num_map_steps(t,time[step],dt);
      break;
    case euler:			// Euler method
      dt = deltat;
      nstep = num_euler_steps(t,time[step],&dt);
      break;
    }

    // loop over individual steps
    for (k = 0; k < nstep; k++) {

      // interpolate the covar functions for the covariates
      table_lookup(&covariate_table,t,cov);

      // loop over replicates
      double *ap, *pm, *xm, *ps = REAL(params);

      for (j = 0, pm = ps, xm = xt; j < nreps; j++, pm += npars, xm += nvars) {

        switch (mode) {

        case Rfun: {

          SEXP ans, nm;

          if (first) {

            PROTECT(ans = eval_call(fn,args,&t,&dt,xm,nvars,pm,npars,cov,ncovars));

            PROTECT(nm = GET_NAMES(ans));
            if (invalid_names(nm))
              errorcall(R_NilValue,"'rprocess' must return a named numeric vector.");
            pidx = INTEGER(PROTECT(matchnames(Snames,nm,"state variables")));

	    nprotect += 3;

            ap = REAL(AS_NUMERIC(ans));
            for (i = 0; i < nvars; i++) xm[pidx[i]] = ap[i];

	    first = 0;

          } else {

            PROTECT(ans = eval_call(fn,args,&t,&dt,xm,nvars,pm,npars,cov,ncovars));

            ap = REAL(AS_NUMERIC(ans));
            for (i = 0; i < nvars; i++) xm[pidx[i]] = ap[i];

            UNPROTECT(1);

          }

        }

          break;

        case native: case regNative: {

          (*ff)(xm,pm,sidx,pidx,cidx,cov,t,dt);

        }

          break;

        default: // # nocov

          errorcall(R_NilValue,"unrecognized 'mode' %d",mode); // # nocov

        }

      }

      t += dt;

      if ((method == euler) && (k == nstep-2)) { // penultimate step
        dt = time[step]-t;
        t = time[step]-dt;
      }

    }

    if (step < ntimes-1)
      memcpy(xt+nvars*nreps,xt,nreps*nvars*sizeof(double));

  }

  // clean up
  switch (mode) {

  case native: case regNative: {

    PutRNGstate();
    unset_pomp_userdata();

  }

    break;

  case Rfun: default:

    break;

  }

  UNPROTECT(nprotect);
  return X;

}

int num_euler_steps (double t1, double t2, double *dt) {
  double tol = sqrt(DOUBLE_EPS);
  int nstep;
  // nstep will be the number of Euler steps to take in going from t1 to t2.
  // note also that the stepsize changes.
  // this choice is meant to be conservative
  // (i.e., so that the actual dt does not exceed the specified dt
  // by more than the relative tolerance 'tol')
  // and to counteract roundoff error.
  // It seems to work well, but is not guaranteed:
  // suggestions would be appreciated.

  if (t1 >= t2) {
    *dt = 0.0;
    nstep = 0;
  } else if (t1+*dt >= t2) {
    *dt = t2-t1;
    nstep = 1;
  } else {
    nstep = (int) ceil((t2-t1)/(*dt)/(1+tol));
    *dt = (t2-t1)/((double) nstep);
  }
  return nstep;
}

int num_map_steps (double t1, double t2, double dt) {
  double tol = sqrt(DOUBLE_EPS);
  int nstep;
  // nstep will be the number of discrete-time steps to take in going from t1 to t2.
  nstep = (int) floor((t2-t1)/dt/(1-tol));
  return (nstep > 0) ? nstep : 0;
}
