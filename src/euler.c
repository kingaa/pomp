// dear emacs, please treat this as -*- C++ -*-

#include "pomp_internal.h"

// take nstep Euler-Poisson steps from t1 to t2
static void euler_simulator (euler_step_sim *estep,
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
  int k, p, step;
  double dt;

  struct lookup_table covariate_table = {covlen, covdim, 0, time_table, covar_table};
  
  GetRNGstate();		// initialize R's pseudorandom number generator

  for (p = 0; p < nrep; p++)
    for (k = 0; k < nvar; k++) 
      x[k+nvar*p] = xstart[k+nvar*p]; // copy the start values into the result array
  
  for (step = 1, t = times[0]; step < ntimes; step++) {	// loop over times

    R_CheckUserInterrupt();	// check for user interrupt

    for (p = 0; p < nrep; p++) {
      xp = &x[nvar*(p+nrep*step)];
      for (k = 0; k < nvar; k++) // copy in the previous values of the state variables
	xp[k] = x[k+nvar*(p+nrep*(step-1))];
      for (k = 0; k < nzero; k++)
	xp[zeroindex[k]] = 0.0;    // set some variables to zero
    }

    dt = *deltat;
    if (t+dt > times[step]) { // if the next step would carry us beyond times[step], reduce dt
      dt = times[step] - t; 
    }

    while (1) {			// loop over Euler steps
      
      // interpolate the covar functions for the covariates
      if (covdim > 0) 
	table_lookup(&covariate_table,t,covar_fn,0);

      for (p = 0; p < nrep; p++) { // loop over replicates
      
	pp = &params[npar*p];
	xp = &x[nvar*(p+nrep*step)];

	(*estep)(xp,pp,stateindex,parindex,covindex,covdim,covar_fn,t,dt);

      }

      t += dt;				    // advance time
      if (t >= times[step]) {		    // finished with Euler steps
	t = times[step];
	break;
      } else if (t > times[step]) { // if the next step would carry us beyond times[step], reduce dt
	dt = times[step] - t; 
      }
    }
  }

  PutRNGstate();		// finished with R's random number generator

}

// these global objects will pass the needed information to the user-defined function (see 'default_euler_step_fn')
// each of these is allocated once, globally, and refilled many times
static SEXP euler_step_Xvec;	// state variable vector
static SEXP euler_step_Pvec;	// parameter vector
static SEXP euler_step_Cvec;	// covariate vector
static SEXP euler_step_time;	// time
static SEXP euler_step_dt;	// stepsize
static int  euler_step_nvar;	// number of state variables
static int  euler_step_npar;	// number of parameters
static SEXP euler_step_envir;	// function's environment
static SEXP euler_step_fcall;	// function call

#define XVEC    (euler_step_Xvec)
#define PVEC    (euler_step_Pvec)
#define CVEC    (euler_step_Cvec)
#define TIME    (euler_step_time)
#define DT      (euler_step_dt)
#define NVAR    (euler_step_nvar)
#define NPAR    (euler_step_npar)
#define RHO     (euler_step_envir)
#define FCALL   (euler_step_fcall)

// this is the euler step function that is evaluated when the user supplies an R function
// (and not a native routine)
// Note that stateindex, parindex, covindex are ignored.
static void default_euler_step_fn (double *x, const double *p, 
				   const int *stateindex, const int *parindex, const int *covindex,
				   int ncovar, const double *covar,
				   double t, double dt)
{
  int nprotect = 0;
  int k;
  double *xp;
  SEXP ans;
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
  xp = REAL(AS_NUMERIC(ans));
  for (k = 0; k < NVAR; k++) x[k] = xp[k];
  UNPROTECT(nprotect);
}


SEXP euler_model_simulator (SEXP func, 
			    SEXP xstart, SEXP times, SEXP params, 
			    SEXP dt, 
			    SEXP statenames, SEXP paramnames, SEXP covarnames, SEXP zeronames,
			    SEXP tcovar, SEXP covar, SEXP args) 
{
  int nprotect = 0;
  int *dim, xdim[3], ndim[7];
  int nvar, npar, nrep, ntimes;
  int covlen, covdim;
  int nstates = LENGTH(statenames);
  int nparams = LENGTH(paramnames);
  int ncovars = LENGTH(covarnames);
  int nzeros = LENGTH(zeronames);
  euler_step_sim *ff = NULL;
  SEXP X, pindex, sindex, cindex, zindex;
  SEXP fn, Xnames, Pnames, Cnames;

  dim = INTEGER(GET_DIM(xstart)); nvar = dim[0]; nrep = dim[1];
  dim = INTEGER(GET_DIM(params)); npar = dim[0];
  dim = INTEGER(GET_DIM(covar)); covlen = dim[0]; covdim = dim[1];
  ntimes = LENGTH(times);

  PROTECT(Xnames = GET_ROWNAMES(GET_DIMNAMES(xstart))); nprotect++;
  PROTECT(Pnames = GET_ROWNAMES(GET_DIMNAMES(params))); nprotect++;
  PROTECT(Cnames = GET_COLNAMES(GET_DIMNAMES(covar))); nprotect++;

  if (inherits(func,"NativeSymbol")) {
    ff = (euler_step_sim *) R_ExternalPtrAddr(func);
  } else if (isFunction(func)) {
    PROTECT(fn = func); nprotect++;
    PROTECT(RHO = (CLOENV(fn))); nprotect++;
    NVAR = nvar;			// for internal use
    NPAR = npar;			// for internal use
    PROTECT(DT = NEW_NUMERIC(1)); nprotect++;	// for internal use
    PROTECT(TIME = NEW_NUMERIC(1)); nprotect++;	// for internal use
    PROTECT(XVEC = NEW_NUMERIC(nvar)); nprotect++; // for internal use
    PROTECT(PVEC = NEW_NUMERIC(npar)); nprotect++; // for internal use
    PROTECT(CVEC = NEW_NUMERIC(covdim)); nprotect++; // for internal use
    SET_NAMES(XVEC,Xnames); // make sure the names attribute is copied
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
    ff = (euler_step_sim *) default_euler_step_fn;
  } else {
    UNPROTECT(nprotect);
    error("illegal input: supplied function must be either an R function or a compiled native function");
  }

  xdim[0] = nvar; xdim[1] = nrep; xdim[2] = ntimes;
  PROTECT(X = makearray(3,xdim)); nprotect++;
  setrownames(X,Xnames,3);
  if (nstates>0) {
    PROTECT(sindex = MATCHROWNAMES(xstart,statenames)); nprotect++;
  } else {
    PROTECT(sindex = NEW_INTEGER(0)); nprotect++;
  }
  if (nparams>0) {
    PROTECT(pindex = MATCHROWNAMES(params,paramnames)); nprotect++;
  } else {
    PROTECT(pindex = NEW_INTEGER(0)); nprotect++;
  }
  if (ncovars>0) {
    PROTECT(cindex = MATCHCOLNAMES(covar,covarnames)); nprotect++;
  } else {
    PROTECT(cindex = NEW_INTEGER(0)); nprotect++;
  }
  if (nzeros>0) {
    PROTECT(zindex = MATCHROWNAMES(xstart,zeronames)); nprotect++;
  } else {
    PROTECT(zindex = NEW_INTEGER(0)); nprotect++;
  }
  ndim[0] = nvar; ndim[1] = npar; ndim[2] = nrep; ndim[3] = ntimes; 
  ndim[4] = covlen; ndim[5] = covdim; ndim[6] = nzeros;
  euler_simulator(ff,REAL(X),REAL(xstart),REAL(times),REAL(params),
		  ndim,REAL(dt),
		  INTEGER(sindex),INTEGER(pindex),INTEGER(cindex),INTEGER(zindex),
		  REAL(tcovar),REAL(covar));
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

// take nstep Euler-Poisson steps from t1 to t2
static void euler_densities (euler_step_pdf *estep,
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

    // interpolate the covar functions for the covariates
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

// these global objects will pass the needed information to the user-defined function (see 'default_euler_dens_fn')
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
static void default_euler_dens_fn (double *f, double *x1, double *x2, double t1, double t2, const double *p, 
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
  euler_step_pdf *ff = NULL;
  SEXP F, pindex, sindex, cindex;
  SEXP fn, Xnames, Pnames, Cnames;

  dim = INTEGER(GET_DIM(x)); nvar = dim[0]; nrep = dim[1];
  dim = INTEGER(GET_DIM(params)); npar = dim[0];
  dim = INTEGER(GET_DIM(covar)); covlen = dim[0]; covdim = dim[1];
  ntimes = LENGTH(times);

  PROTECT(Xnames = GET_ROWNAMES(GET_DIMNAMES(x))); nprotect++;
  PROTECT(Pnames = GET_ROWNAMES(GET_DIMNAMES(params))); nprotect++;
  PROTECT(Cnames = GET_COLNAMES(GET_DIMNAMES(covar))); nprotect++;

  if (inherits(func,"NativeSymbol")) {
    ff = (euler_step_pdf *) R_ExternalPtrAddr(func);
  } else if (isFunction(func)) {
    PROTECT(fn = func); nprotect++;
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
    ff = (euler_step_pdf *) default_euler_dens_fn;
  } else {
    UNPROTECT(nprotect);
    error("illegal input: supplied function must be either an R function or a compiled native function");
  }

  fdim[0] = nrep; fdim[1] = ntimes-1;
  PROTECT(F = makearray(2,fdim)); nprotect++;
  if (nstates>0) {
    PROTECT(sindex = MATCHROWNAMES(x,statenames)); nprotect++;
  } else {
    PROTECT(sindex = NEW_INTEGER(0)); nprotect++;
  }
  if (nparams>0) {
    PROTECT(pindex = MATCHROWNAMES(params,paramnames)); nprotect++;
  } else {
    PROTECT(pindex = NEW_INTEGER(0)); nprotect++;
  }
  if (ncovars>0) {
    PROTECT(cindex = MATCHCOLNAMES(covar,covarnames)); nprotect++;
  } else {
    PROTECT(cindex = NEW_INTEGER(0)); nprotect++;
  }
  ndim[0] = nvar; ndim[1] = npar; ndim[2] = nrep; ndim[3] = ntimes; 
  ndim[4] = covlen; ndim[5] = covdim;
  euler_densities(ff,REAL(F),REAL(x),REAL(times),REAL(params),
		  ndim,INTEGER(sindex),INTEGER(pindex),INTEGER(cindex),
		  REAL(tcovar),REAL(covar),INTEGER(log));
  UNPROTECT(nprotect);
  return F;
}

#undef X1VEC
#undef X2VEC
#undef PVEC
#undef CVEC
#undef TIME
#undef NVAR
#undef NPAR
#undef RHO
#undef FCALL
