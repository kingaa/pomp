// dear emacs, please treat this as -*- C++ -*-

#include "pomp_internal.h"
#include <R_ext/Constants.h>

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

// take Euler-Poisson steps of size at most deltat from t1 to t2
static void euler_simulator (pomp_onestep_sim *estep,
			     double *x, double *xstart, double *times, double *params, 
			     int *ndim, double *deltat,
			     int *stateindex, int *parindex, int *covindex, int *zeroindex,
			     double *time_table, double *covar_table)
{
  double t, *xp, *pp;
  int nvar = ndim[0];
  int npar = ndim[1];
  int nrep = ndim[2];
  int ntimes = ndim[3];
  int covlen = ndim[4];
  int covdim = ndim[5];
  int nzero = ndim[6];
  double covar_fn[covdim];
  int j, k, p, step, nstep = 0;
  double dt;

  struct lookup_table covariate_table = {covlen, covdim, 0, time_table, covar_table};

  // copy the start values into the result array
  for (p = 0; p < nrep; p++)
    for (k = 0; k < nvar; k++) 
      x[k+nvar*p] = xstart[k+nvar*p];
  
  // loop over times
  for (step = 1; step < ntimes; step++) {

    R_CheckUserInterrupt();

    t = times[step-1];
    dt = *deltat;

    if (t > times[step]) {
      error("'times' is not an increasing sequence");
    }

    nstep = num_euler_steps(t,times[step],&dt);

    for (p = 0; p < nrep; p++) {
      xp = &x[nvar*(p+nrep*step)];
      // copy in the previous values of the state variables
      for (k = 0; k < nvar; k++)
	xp[k] = x[k+nvar*(p+nrep*(step-1))];
      // set some variables to zero
      for (k = 0; k < nzero; k++)
	xp[zeroindex[k]] = 0.0;
    }

    for (j = 0; j < nstep; j++) { // loop over Euler steps

      // interpolate the covar functions for the covariates
      if (covdim > 0) 
	table_lookup(&covariate_table,t,covar_fn,0);

      for (p = 0; p < nrep; p++) { // loop over replicates
      
	pp = &params[npar*p];
	xp = &x[nvar*(p+nrep*step)];

	(*estep)(xp,pp,stateindex,parindex,covindex,covdim,covar_fn,t,dt);

      }

      t += dt;

      if (j == nstep-2) {	// penultimate step
	dt = times[step]-t;
	t = times[step]-dt;
      }

    }
  }
}

// take discrete steps of size deltat from t1 to t2
static void discrete_time_simulator (pomp_onestep_sim *estep,
				     double *x, double *xstart, double *times, double *params, 
				     int *ndim, double *deltat,
				     int *stateindex, int *parindex, int *covindex, int *zeroindex,
				     double *time_table, double *covar_table)
{
  double t, *xp, *pp;
  int nvar = ndim[0];
  int npar = ndim[1];
  int nrep = ndim[2];
  int ntimes = ndim[3];
  int covlen = ndim[4];
  int covdim = ndim[5];
  int nzero = ndim[6];
  double covar_fn[covdim];
  int j, k, p, step, nstep = 0;
  double dt;

  struct lookup_table covariate_table = {covlen, covdim, 0, time_table, covar_table};

  // copy the start values into the result array
  for (p = 0; p < nrep; p++)
    for (k = 0; k < nvar; k++) 
      x[k+nvar*p] = xstart[k+nvar*p];
  
  dt = *deltat;
  t = times[0];

  // loop over times
  for (step = 1; step < ntimes; step++) {

    R_CheckUserInterrupt();

    nstep = num_map_steps(t,times[step],dt);

    for (p = 0; p < nrep; p++) {
      xp = &x[nvar*(p+nrep*step)];
      // copy in the previous values of the state variables
      for (k = 0; k < nvar; k++)
	xp[k] = x[k+nvar*(p+nrep*(step-1))];
      // set some variables to zero
      for (k = 0; k < nzero; k++)
	xp[zeroindex[k]] = 0.0;
    }

    for (j = 0; j < nstep; j++) { // loop over steps

      // interpolate the covar functions for the covariates
      if (covdim > 0) 
	table_lookup(&covariate_table,t,covar_fn,0);

      for (p = 0; p < nrep; p++) { // loop over replicates
      
	pp = &params[npar*p];
	xp = &x[nvar*(p+nrep*step)];

	(*estep)(xp,pp,stateindex,parindex,covindex,covdim,covar_fn,t,dt);

      }

      t += dt;

    }
  }
}

// take one step from t1 to t2
static void onestep_simulator (pomp_onestep_sim *estep,
			       double *x, double *xstart, double *times, double *params, 
			       int *ndim, 
			       int *stateindex, int *parindex, int *covindex, int *zeroindex,
			       double *time_table, double *covar_table)
{
  double t, *xp, *pp;
  int nvar = ndim[0];
  int npar = ndim[1];
  int nrep = ndim[2];
  int ntimes = ndim[3];
  int covlen = ndim[4];
  int covdim = ndim[5];
  int nzero = ndim[6];
  double covar_fn[covdim];
  int k, p, step;
  double dt;

  struct lookup_table covariate_table = {covlen, covdim, 0, time_table, covar_table};

  // copy the start values into the result array
  for (p = 0; p < nrep; p++)
    for (k = 0; k < nvar; k++) 
      x[k+nvar*p] = xstart[k+nvar*p];
  
  // loop over times
  for (step = 1; step < ntimes; step++) {

    R_CheckUserInterrupt();

    t = times[step-1];
    dt = times[step]-t;

    // interpolate the covar functions for the covariates
    if (covdim > 0) 
      table_lookup(&covariate_table,t,covar_fn,0);
    
    for (p = 0; p < nrep; p++) {
      xp = &x[nvar*(p+nrep*step)];
      // copy in the previous values of the state variables
      for (k = 0; k < nvar; k++)
	xp[k] = x[k+nvar*(p+nrep*(step-1))];
      // set some variables to zero
      for (k = 0; k < nzero; k++)
	xp[zeroindex[k]] = 0.0;

      if (t < times[step]) { //  call the simulator only if time elapses
	pp = &params[npar*p];
	xp = &x[nvar*(p+nrep*step)];
	(*estep)(xp,pp,stateindex,parindex,covindex,covdim,covar_fn,t,dt);
      }

    }
  }
}

// these global objects will pass the needed information to the user-defined function (see 'default_onestep_sim_fn')
// each of these is allocated once, globally, and refilled many times
static SEXP _onestep_internal_Xvec;	// state variable vector
static SEXP _onestep_internal_Pvec;	// parameter vector
static SEXP _onestep_internal_Cvec;	// covariate vector
static SEXP _onestep_internal_time;	// time
static SEXP _onestep_internal_dt;	// stepsize
static int  _onestep_internal_nvar;	// number of state variables
static int  _onestep_internal_npar;	// number of parameters
static SEXP _onestep_internal_envir;	// function's environment
static SEXP _onestep_internal_fcall;	// function call
static int  _onestep_internal_first;	// first evaluation?
static SEXP _onestep_internal_vnames;	// names of state variables
static int *_onestep_internal_vindex;	// indices of state variables

#define FIRST   (_onestep_internal_first)
#define VNAMES  (_onestep_internal_vnames)
#define VINDEX  (_onestep_internal_vindex)
#define XVEC    (_onestep_internal_Xvec)
#define PVEC    (_onestep_internal_Pvec)
#define CVEC    (_onestep_internal_Cvec)
#define TIME    (_onestep_internal_time)
#define DT      (_onestep_internal_dt)
#define NVAR    (_onestep_internal_nvar)
#define NPAR    (_onestep_internal_npar)
#define RHO     (_onestep_internal_envir)
#define FCALL   (_onestep_internal_fcall)

// this is the euler step function that is evaluated when the user supplies an R function
// (and not a native routine)
// Note that stateindex, parindex, covindex are ignored.
static void default_onestep_sim_fn (double *x, const double *p, 
				    const int *stateindex, const int *parindex, const int *covindex,
				    int ncovar, const double *covar,
				    double t, double dt)
{
  int nprotect = 0;
  int *op, k;
  double *xp;
  SEXP ans, nm, idx;
  xp = REAL(XVEC);
  for (k = 0; k < NVAR; k++) xp[k] = x[k];
  xp = REAL(PVEC);
  for (k = 0; k < NPAR; k++) xp[k] = p[k];
  xp = REAL(CVEC);
  for (k = 0; k < ncovar; k++) xp[k] = covar[k];
  xp = REAL(TIME);
  xp[0] = t;
  xp = REAL(DT);
  xp[0] = dt;

  PROTECT(ans = eval(FCALL,RHO)); nprotect++; // evaluate the call

  if (FIRST) {
    if (LENGTH(ans) != NVAR) {
      UNPROTECT(nprotect);
      error("user 'step.fun' returns a vector of %d states but %d are expected: compare initial conditions?",
	    LENGTH(ans),NVAR);
    }
    PROTECT(nm = GET_NAMES(ans)); nprotect++;
    if (!isNull(nm)) {	   // match names against names from data slot
      PROTECT(idx = matchnames(VNAMES,nm)); nprotect++;
      op = INTEGER(idx);
      for (k = 0; k < NVAR; k++) VINDEX[k] = op[k];
    } else {
      VINDEX = 0;
    }
    FIRST = 0;
  }

  xp = REAL(AS_NUMERIC(ans));
  if (VINDEX != 0) {
    for (k = 0; k < NVAR; k++) x[VINDEX[k]] = xp[k];
  } else {
    for (k = 0; k < NVAR; k++) x[k] = xp[k];
  }

  UNPROTECT(nprotect);
}


SEXP euler_model_simulator (SEXP func, 
			    SEXP xstart, SEXP times, SEXP params, 
			    SEXP dt, SEXP method,
			    SEXP statenames, SEXP paramnames, SEXP covarnames, SEXP zeronames,
			    SEXP tcovar, SEXP covar, SEXP args) 
{
  int nprotect = 0;
  int use_native = 0;
  int *dim, xdim[3], ndim[7];
  int nvar, npar, nrep, ntimes;
  int covlen, covdim;
  int nstates = LENGTH(statenames);
  int nparams = LENGTH(paramnames);
  int ncovars = LENGTH(covarnames);
  int nzeros = LENGTH(zeronames);
  pomp_onestep_sim *ff = NULL;
  SEXP X, pindex, sindex, cindex, zindex;
  int *sidx, *pidx, *cidx, *zidx;
  SEXP fn, Pnames, Cnames;
  int meth = 0;

  dim = INTEGER(GET_DIM(xstart)); nvar = dim[0]; nrep = dim[1];
  dim = INTEGER(GET_DIM(params)); npar = dim[0];
  dim = INTEGER(GET_DIM(covar)); covlen = dim[0]; covdim = dim[1];
  ntimes = LENGTH(times);

  meth = *(INTEGER(AS_INTEGER(method))); // 0 = Euler, 1 = one-step, 2 = fixed step

  PROTECT(VNAMES = GET_ROWNAMES(GET_DIMNAMES(xstart))); nprotect++;
  PROTECT(Pnames = GET_ROWNAMES(GET_DIMNAMES(params))); nprotect++;
  PROTECT(Cnames = GET_COLNAMES(GET_DIMNAMES(covar))); nprotect++;

  PROTECT(fn = pomp_fun_handler(func,&use_native)); nprotect++;

  if (use_native) {
    ff = (pomp_onestep_sim *) R_ExternalPtrAddr(fn);
    VINDEX = 0;
  } else {
    PROTECT(RHO = (CLOENV(fn))); nprotect++;
    NVAR = nvar;			// for internal use
    NPAR = npar;			// for internal use
    PROTECT(DT = NEW_NUMERIC(1)); nprotect++;	// for internal use
    PROTECT(TIME = NEW_NUMERIC(1)); nprotect++;	// for internal use
    PROTECT(XVEC = NEW_NUMERIC(nvar)); nprotect++; // for internal use
    PROTECT(PVEC = NEW_NUMERIC(npar)); nprotect++; // for internal use
    PROTECT(CVEC = NEW_NUMERIC(covdim)); nprotect++; // for internal use
    SET_NAMES(XVEC,VNAMES); // make sure the names attribute is copied
    SET_NAMES(PVEC,Pnames); // make sure the names attribute is copied
    SET_NAMES(CVEC,Cnames); // make sure the names attribute is copied
    // set up the function call
    PROTECT(FCALL = LCONS(CVEC,args)); nprotect++;
    SET_TAG(FCALL,install("covars"));
    PROTECT(FCALL = LCONS(DT,FCALL)); nprotect++;
    SET_TAG(FCALL,install("delta.t"));
    PROTECT(FCALL = LCONS(PVEC,FCALL)); nprotect++;
    SET_TAG(FCALL,install("params"));
    PROTECT(FCALL = LCONS(TIME,FCALL)); nprotect++;
    SET_TAG(FCALL,install("t"));
    PROTECT(FCALL = LCONS(XVEC,FCALL)); nprotect++;
    SET_TAG(FCALL,install("x"));
    PROTECT(FCALL = LCONS(fn,FCALL)); nprotect++;
    ff = (pomp_onestep_sim *) default_onestep_sim_fn;
    VINDEX = (int *) R_alloc(nvar,sizeof(int));
    FIRST = 1;
  }

  xdim[0] = nvar; xdim[1] = nrep; xdim[2] = ntimes;
  PROTECT(X = makearray(3,xdim)); nprotect++;
  setrownames(X,VNAMES,3);

  if (nstates>0) {
    PROTECT(sindex = MATCHROWNAMES(xstart,statenames)); nprotect++;
    sidx = INTEGER(sindex);
  } else {
    sidx = 0;
  }
  if (nparams>0) {
    PROTECT(pindex = MATCHROWNAMES(params,paramnames)); nprotect++;
    pidx = INTEGER(pindex);
  } else {
    pidx = 0;
  }
  if (ncovars>0) {
    PROTECT(cindex = MATCHCOLNAMES(covar,covarnames)); nprotect++;
    cidx = INTEGER(cindex);
  } else {
    cidx = 0;
  }
  if (nzeros>0) {
    PROTECT(zindex = MATCHROWNAMES(xstart,zeronames)); nprotect++;
    zidx = INTEGER(zindex);
  } else {
    zidx = 0;
  }

  ndim[0] = nvar; ndim[1] = npar; ndim[2] = nrep; ndim[3] = ntimes; 
  ndim[4] = covlen; ndim[5] = covdim; ndim[6] = nzeros;

  if (use_native) GetRNGstate();

  switch (meth) {
  case 0:
    euler_simulator(ff,REAL(X),REAL(xstart),REAL(times),REAL(params),
		    ndim,REAL(dt),sidx,pidx,cidx,zidx,
		    REAL(tcovar),REAL(covar));
    break;
  case 1:
    onestep_simulator(ff,REAL(X),REAL(xstart),REAL(times),REAL(params),
		      ndim,sidx,pidx,cidx,zidx,
		      REAL(tcovar),REAL(covar));
    break;
  case 2:
    discrete_time_simulator(ff,REAL(X),REAL(xstart),REAL(times),REAL(params),
			    ndim,REAL(dt),sidx,pidx,cidx,zidx,
			    REAL(tcovar),REAL(covar));
    break;
  }
  
  if (use_native) PutRNGstate();

  VINDEX = 0;

  UNPROTECT(nprotect);
  return X;
}

#undef XVEC
#undef PVEC
#undef CVEC
#undef TIME
#undef DT
#undef NVAR
#undef NPAR
#undef RHO
#undef FCALL
#undef FIRST
#undef VNAMES
#undef VINDEX

// compute pdf of a sequence of Euler steps
static void euler_densities (pomp_onestep_pdf *estep,
			     double *f, 
			     double *x, double *times, double *params, 
			     int *ndim,
			     int *stateindex, int *parindex, int *covindex,
			     double *time_table, double *covar_table,
			     int *give_log)
{
  double *x1p, *x2p, *pp, *fp, t1, t2;
  int nvar = ndim[0];
  int npar = ndim[1];
  int nrep = ndim[2];
  int ntimes = ndim[3];
  int covlen = ndim[4];
  int covdim = ndim[5];
  double covar_fn[covdim];
  int p, step;

  // set up the covariate table
  struct lookup_table covariate_table = {covlen, covdim, 0, time_table, covar_table};
  
  for (step = 0; step < ntimes-1; step++) { // loop over times

    R_CheckUserInterrupt();	// check for user interrupt

    t1 = times[step];
    t2 = times[step+1];

    // interpolate the covariates at time t1
    if (covdim > 0) 
      table_lookup(&covariate_table,t1,covar_fn,0);
    
    for (p = 0; p < nrep; p++) { // loop over replicates
      
      fp = &f[p+nrep*step];
      x1p = &x[nvar*(p+nrep*step)];
      x2p = &x[nvar*(p+nrep*(step+1))];
      pp = &params[npar*p];
      
      (*estep)(fp,x1p,x2p,t1,t2,pp,
	       stateindex,parindex,covindex,
	       covdim,covar_fn);
      
      if (!(*give_log)) *fp = exp(*fp);
      
    }
  }
}

// these global objects will pass the needed information to the user-defined function (see 'default_onestep_dens_fn')
// each of these is allocated once, globally, and refilled many times
static SEXP euler_dens_Xvec1;	// state variable vector
static SEXP euler_dens_Xvec2;	// state variable vector
static SEXP euler_dens_Pvec;	// parameter vector
static SEXP euler_dens_Cvec;	// covariate vector
static SEXP euler_dens_time1;	// time1
static SEXP euler_dens_time2;	// time2
static int  euler_dens_nvar;	// number of state variables
static int  euler_dens_npar;	// number of parameters
static SEXP euler_dens_envir;	// function's environment
static SEXP euler_dens_fcall;	// function call

#define X1VEC   (euler_dens_Xvec1)
#define X2VEC   (euler_dens_Xvec2)
#define PVEC    (euler_dens_Pvec)
#define CVEC    (euler_dens_Cvec)
#define TIME1   (euler_dens_time1)
#define TIME2   (euler_dens_time2)
#define NVAR    (euler_dens_nvar)
#define NPAR    (euler_dens_npar)
#define RHO     (euler_dens_envir)
#define FCALL   (euler_dens_fcall)

// this is the euler dens function that is evaluated when the user supplies an R function
// (and not a native routine)
// Note that stateindex, parindex, covindex are ignored.
static void default_onestep_dens_fn (double *f, double *x1, double *x2, double t1, double t2, const double *p, 
				     const int *stateindex, const int *parindex, const int *covindex,
				     int ncovar, const double *covar)
{
  int nprotect = 0;
  int k;
  double *xp;
  SEXP ans;
  xp = REAL(X1VEC);
  for (k = 0; k < NVAR; k++) xp[k] = x1[k];
  xp = REAL(X2VEC);
  for (k = 0; k < NVAR; k++) xp[k] = x2[k];
  xp = REAL(PVEC);
  for (k = 0; k < NPAR; k++) xp[k] = p[k];
  xp = REAL(CVEC);
  for (k = 0; k < ncovar; k++) xp[k] = covar[k];
  xp = REAL(TIME1); xp[0] = t1;
  xp = REAL(TIME2); xp[0] = t2;

  PROTECT(ans = eval(FCALL,RHO)); nprotect++; // evaluate the call

  xp = REAL(AS_NUMERIC(ans));
  f[0] = xp[0];
  UNPROTECT(nprotect);
}


SEXP euler_model_density (SEXP func, 
			  SEXP x, SEXP times, SEXP params, 
			  SEXP statenames, SEXP paramnames, SEXP covarnames,
			  SEXP tcovar, SEXP covar, SEXP log, SEXP args) 
{
  int nprotect = 0;
  int *dim, fdim[2], ndim[6];
  int nvar, npar, nrep, ntimes;
  int covlen, covdim;
  int nstates = LENGTH(statenames);
  int nparams = LENGTH(paramnames);
  int ncovars = LENGTH(covarnames);
  pomp_onestep_pdf *ff = NULL;
  SEXP F, pindex, sindex, cindex;
  int *pidx, *sidx, *cidx;
  SEXP fn, Xnames, Pnames, Cnames;
  int use_native;

  dim = INTEGER(GET_DIM(x)); nvar = dim[0]; nrep = dim[1];
  dim = INTEGER(GET_DIM(params)); npar = dim[0];
  dim = INTEGER(GET_DIM(covar)); covlen = dim[0]; covdim = dim[1];
  ntimes = LENGTH(times);

  PROTECT(Xnames = GET_ROWNAMES(GET_DIMNAMES(x))); nprotect++;
  PROTECT(Pnames = GET_ROWNAMES(GET_DIMNAMES(params))); nprotect++;
  PROTECT(Cnames = GET_COLNAMES(GET_DIMNAMES(covar))); nprotect++;

  PROTECT(fn = pomp_fun_handler(func,&use_native)); nprotect++;

  if (use_native) {
    ff = (pomp_onestep_pdf *) R_ExternalPtrAddr(fn);
  } else {
    PROTECT(RHO = (CLOENV(fn))); nprotect++;
    NVAR = nvar;			// for internal use
    NPAR = npar;			// for internal use
    PROTECT(TIME1 = NEW_NUMERIC(1)); nprotect++;	// for internal use
    PROTECT(TIME2 = NEW_NUMERIC(1)); nprotect++;	// for internal use
    PROTECT(X1VEC = NEW_NUMERIC(nvar)); nprotect++; // for internal use
    PROTECT(X2VEC = NEW_NUMERIC(nvar)); nprotect++; // for internal use
    PROTECT(PVEC = NEW_NUMERIC(npar)); nprotect++; // for internal use
    PROTECT(CVEC = NEW_NUMERIC(covdim)); nprotect++; // for internal use
    SET_NAMES(X1VEC,Xnames); // make sure the names attribute is copied
    SET_NAMES(X2VEC,Xnames); // make sure the names attribute is copied
    SET_NAMES(PVEC,Pnames); // make sure the names attribute is copied
    SET_NAMES(CVEC,Cnames); // make sure the names attribute is copied
    // set up the function call
    PROTECT(FCALL = LCONS(CVEC,args)); nprotect++;
    SET_TAG(FCALL,install("covars"));
    PROTECT(FCALL = LCONS(PVEC,FCALL)); nprotect++;
    SET_TAG(FCALL,install("params"));
    PROTECT(FCALL = LCONS(TIME2,FCALL)); nprotect++;
    SET_TAG(FCALL,install("t2"));
    PROTECT(FCALL = LCONS(TIME1,FCALL)); nprotect++;
    SET_TAG(FCALL,install("t1"));
    PROTECT(FCALL = LCONS(X2VEC,FCALL)); nprotect++;
    SET_TAG(FCALL,install("x2"));
    PROTECT(FCALL = LCONS(X1VEC,FCALL)); nprotect++;
    SET_TAG(FCALL,install("x1"));
    PROTECT(FCALL = LCONS(fn,FCALL)); nprotect++;
    ff = (pomp_onestep_pdf *) default_onestep_dens_fn;
  }

  fdim[0] = nrep; fdim[1] = ntimes-1;
  PROTECT(F = makearray(2,fdim)); nprotect++;

  if (nstates>0) {
    PROTECT(sindex = MATCHROWNAMES(x,statenames)); nprotect++;
    sidx = INTEGER(sindex);
  } else {
    sidx = 0;
  }
  if (nparams>0) {
    PROTECT(pindex = MATCHROWNAMES(params,paramnames)); nprotect++;
    pidx = INTEGER(pindex);
  } else {
    pidx = 0;
  }
  if (ncovars>0) {
    PROTECT(cindex = MATCHCOLNAMES(covar,covarnames)); nprotect++;
    cidx = INTEGER(cindex);
  } else {
    cidx = 0;
  }

  ndim[0] = nvar; ndim[1] = npar; ndim[2] = nrep; ndim[3] = ntimes; 
  ndim[4] = covlen; ndim[5] = covdim;

  euler_densities(ff,REAL(F),REAL(x),REAL(times),REAL(params),
		  ndim,sidx,pidx,cidx,
		  REAL(tcovar),REAL(covar),INTEGER(log));

  UNPROTECT(nprotect);
  return F;
}

#undef X1VEC
#undef X2VEC
#undef PVEC
#undef CVEC
#undef TIME1
#undef TIME2
#undef NVAR
#undef NPAR
#undef RHO
#undef FCALL
