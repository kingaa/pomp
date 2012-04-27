// dear emacs, please treat this as -*- C++ -*-

#include "pomp_internal.h"
#include <R_ext/Constants.h>

SEXP euler_model_simulator (SEXP func, 
                            SEXP xstart, SEXP times, SEXP params, 
                            SEXP deltat, SEXP method,
                            SEXP statenames, SEXP paramnames, SEXP covarnames, SEXP zeronames,
                            SEXP tcovar, SEXP covar, SEXP args) 
{
  int nprotect = 0;
  int first = 1;
  int mode = -1;
  int meth = 0;
  int use_names;
  int *dim, *posn, xdim[3];
  int nstep, nvars, npars, nreps, ntimes, nzeros, ncovars, covlen;
  pomp_onestep_sim *ff = NULL;
  SEXP X;
  SEXP fn, fcall, ans, rho, nm;
  SEXP Snames, Pnames, Cnames;
  SEXP tvec, xvec, pvec, cvec, dtvec;
  int *sidx, *pidx, *cidx, *zidx;
  double *tp, *xp, *pp, *cp, *dtp, *xt, *ps;
  double *time;
  double *Xt;
  double t, dt;
  int i, j, k, step;

  dim = INTEGER(GET_DIM(xstart)); nvars = dim[0]; nreps = dim[1];
  dim = INTEGER(GET_DIM(params)); npars = dim[0];
  dim = INTEGER(GET_DIM(covar)); covlen = dim[0]; ncovars = dim[1];
  ntimes = LENGTH(times);

  PROTECT(Snames = GET_ROWNAMES(GET_DIMNAMES(xstart))); nprotect++;
  PROTECT(Pnames = GET_ROWNAMES(GET_DIMNAMES(params))); nprotect++;
  PROTECT(Cnames = GET_COLNAMES(GET_DIMNAMES(covar))); nprotect++;

  // set up the covariate table
  struct lookup_table covariate_table = {covlen, ncovars, 0, REAL(tcovar), REAL(covar)};

  // vector for interpolated covariates
  PROTECT(cvec = NEW_NUMERIC(ncovars)); nprotect++;
  SET_NAMES(cvec,Cnames);
  cp = REAL(cvec);

  // indices of accumulator variables
  nzeros = LENGTH(zeronames);
  zidx = INTEGER(PROTECT(matchnames(Snames,zeronames))); nprotect++;

  // extract user function
  PROTECT(fn = pomp_fun_handler(func,&mode)); nprotect++;
  
  // set up
  switch (mode) {

  case 1:			// native code

    // construct state, parameter, covariate, observable indices
    sidx = INTEGER(PROTECT(matchnames(Snames,statenames))); nprotect++;
    pidx = INTEGER(PROTECT(matchnames(Pnames,paramnames))); nprotect++;
    cidx = INTEGER(PROTECT(matchnames(Cnames,covarnames))); nprotect++;

    ff = (pomp_onestep_sim *) R_ExternalPtrAddr(fn);

    break;

  case 0:			// R function

    // get function's environment
    PROTECT(rho = (CLOENV(fn))); nprotect++;

    PROTECT(dtvec = NEW_NUMERIC(1)); nprotect++;
    dtp = REAL(dtvec);
    PROTECT(tvec = NEW_NUMERIC(1)); nprotect++;
    tp = REAL(tvec);
    PROTECT(xvec = NEW_NUMERIC(nvars)); nprotect++;
    SET_NAMES(xvec,Snames);
    xp = REAL(xvec);
    PROTECT(pvec = NEW_NUMERIC(npars)); nprotect++;
    SET_NAMES(pvec,Pnames);
    pp = REAL(pvec);

    // set up the function call
    PROTECT(fcall = LCONS(cvec,args)); nprotect++;
    SET_TAG(fcall,install("covars"));
    PROTECT(fcall = LCONS(dtvec,fcall)); nprotect++;
    SET_TAG(fcall,install("delta.t"));
    PROTECT(fcall = LCONS(pvec,fcall)); nprotect++;
    SET_TAG(fcall,install("params"));
    PROTECT(fcall = LCONS(tvec,fcall)); nprotect++;
    SET_TAG(fcall,install("t"));
    PROTECT(fcall = LCONS(xvec,fcall)); nprotect++;
    SET_TAG(fcall,install("x"));
    PROTECT(fcall = LCONS(fn,fcall)); nprotect++;

    // to hold indices of variables that must be rearranged
    posn = (int *) R_alloc(nvars,sizeof(int));    

    // this is the euler step function that is evaluated when the user supplies an R function
    // (and not a native routine)
    // Note that stateindex, parindex, covindex are ignored.
    void R_step_fn (double *x, const double *p, 
		    const int *stateindex, const int *parindex, const int *covindex,
		    int ncovar, const double *covar,
		    double t, double dt)
    {
      int nprotect = 0;
      int *op, i;
      double *xs;
      SEXP ans, nm;

      for (i = 0; i < nvars; i++) xp[i] = x[i];
      for (i = 0; i < npars; i++) pp[i] = p[i];
      *tp = t;
      *dtp = dt;
      
      if (first) {

	PROTECT(ans = eval(fcall,rho)); nprotect++; // evaluate the call
      
	if (LENGTH(ans) != nvars) {
	  error("user 'step.fun' returns a vector of %d states but %d are expected: compare initial conditions?",
		LENGTH(ans),nvars);
	}

	PROTECT(nm = GET_NAMES(ans)); nprotect++;
	use_names = !isNull(nm);
	if (use_names) {
	  op = INTEGER(PROTECT(matchnames(Snames,nm))); nprotect++;
	  for (i = 0; i < nvars; i++) posn[i] = op[i];
	}

	xs = REAL(AS_NUMERIC(ans));

	first = 0;

      } else {
      
	xs = REAL(AS_NUMERIC(eval(fcall,rho)));

      }

      if (use_names) {
	for (i = 0; i < nvars; i++) x[posn[i]] = xs[i];
      } else {
	for (i = 0; i < nvars; i++) x[i] = xs[i];
      }
      
      UNPROTECT(nprotect);
    }

    sidx = 0;
    pidx = 0;
    cidx = 0;

    ff = (pomp_onestep_sim *) R_step_fn;

    break;

  default:
    error("unrecognized 'mode' in 'euler_simulator'");
    break;
  }

  if (mode==1) {
    set_pomp_userdata(args);
    GetRNGstate();
  }

  // create array to hold results
  xdim[0] = nvars; xdim[1] = nreps; xdim[2] = ntimes;
  PROTECT(X = makearray(3,xdim)); nprotect++;
  setrownames(X,Snames,3);
  Xt = REAL(X);

  // copy the start values into the result array
  xt = REAL(xstart);
  for (j = 0; j < nreps; j++)
    for (i = 0; i < nvars; i++) 
      Xt[i+nvars*j] = xt[i+nvars*j];

  meth = *(INTEGER(AS_INTEGER(method))); // 0 = Euler, 1 = one-step, 2 = fixed step

  // now do computations
  // loop over times
  time = REAL(times);
  t = time[0];

  for (step = 1; step < ntimes; step++) {

    R_CheckUserInterrupt();

    if (t > time[step]) {
      error("'times' is not an increasing sequence");
    }

    switch (meth) {
    case 0:			// Euler method
      dt = *(REAL(deltat));
      nstep = num_euler_steps(t,time[step],&dt);
      break;
    case 1:			// one step 
      dt = time[step]-t;
      nstep = (dt > 0) ? 1 : 0;
      break;
    case 2:			// fixed step
      dt = *(REAL(deltat));
      nstep = num_map_steps(t,time[step],dt);
      break;
    default:
      error("unrecognized 'method' in 'stepwise_simulator'");
    }

    for (j = 0; j < nreps; j++) {
      xt = &Xt[nvars*(j+nreps*step)];
      // copy in the previous values of the state variables
      for (i = 0; i < nvars; i++) xt[i] = Xt[i+nvars*(j+nreps*(step-1))];
      // set some variables to zero 
      for (i = 0; i < nzeros; i++) xt[zidx[i]] = 0.0;
    }

    for (k = 0; k < nstep; k++) { // loop over Euler steps

      // interpolate the covar functions for the covariates
      table_lookup(&covariate_table,t,cp,0);

      for (j = 0, ps = REAL(params); j < nreps; j++, ps += npars) { // loop over replicates
      
	xt = &Xt[nvars*(j+nreps*step)];

	(*ff)(xt,ps,sidx,pidx,cidx,ncovars,cp,t,dt);

      }

      t += dt;

      if ((method == 0) && (k == nstep-2)) { // penultimate step
	dt = time[step]-t;
	t = time[step]-dt;
      }

    }
  }

  if (mode==1) {
    PutRNGstate();
    unset_pomp_userdata();
  }

  UNPROTECT(nprotect);
  return X;
}

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

  if (use_native) set_pomp_userdata(args);

  euler_densities(ff,REAL(F),REAL(x),REAL(times),REAL(params),
		  ndim,sidx,pidx,cidx,
		  REAL(tcovar),REAL(covar),INTEGER(log));

  if (use_native) unset_pomp_userdata();

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

