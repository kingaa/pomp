// dear emacs, please treat this as -*- C++ -*-

#include "pomp_internal.h"
#include <R_ext/Constants.h>

SEXP euler_model_simulator (SEXP func, 
                            SEXP xstart, SEXP times, SEXP params, 
                            SEXP deltat, SEXP method,
                            SEXP statenames, SEXP paramnames, SEXP covarnames, SEXP zeronames,
                            SEXP tcovar, SEXP covar, SEXP args, SEXP gnsi) 
{
  int nprotect = 0;
  int mode = -1;
  int nstep, nvars, npars, nreps, ntimes, nzeros, ncovars, covlen;
  SEXP X;
  SEXP fn, fcall, rho, ans, nm;
  SEXP Snames, Pnames, Cnames;
  SEXP tvec, xvec, pvec, cvec, dtvec;
  int *pidx = 0, *sidx = 0, *cidx = 0, *zidx = 0;
  pomp_onestep_sim *ff = NULL;
  int meth = *(INTEGER(AS_INTEGER(method))); // 0 = Euler, 1 = one-step, 2 = fixed step

  {
    int *dim;
    dim = INTEGER(GET_DIM(xstart)); nvars = dim[0]; nreps = dim[1];
    dim = INTEGER(GET_DIM(params)); npars = dim[0];
    dim = INTEGER(GET_DIM(covar)); covlen = dim[0]; ncovars = dim[1];
    ntimes = LENGTH(times);
  }

  PROTECT(Snames = GET_ROWNAMES(GET_DIMNAMES(xstart))); nprotect++;
  PROTECT(Pnames = GET_ROWNAMES(GET_DIMNAMES(params))); nprotect++;
  PROTECT(Cnames = GET_COLNAMES(GET_DIMNAMES(covar))); nprotect++;

  // set up the covariate table
  struct lookup_table covariate_table = {covlen, ncovars, 0, REAL(tcovar), REAL(covar)};

  // vector for interpolated covariates
  PROTECT(cvec = NEW_NUMERIC(ncovars)); nprotect++;
  SET_NAMES(cvec,Cnames);

  // indices of accumulator variables
  nzeros = LENGTH(zeronames);
  zidx = INTEGER(PROTECT(matchnames(Snames,zeronames))); nprotect++;

  // extract user function
  PROTECT(fn = pomp_fun_handler(func,gnsi,&mode)); nprotect++;
  
  // set up
  switch (mode) {

  case 0:			// R function

    PROTECT(dtvec = NEW_NUMERIC(1)); nprotect++;
    PROTECT(tvec = NEW_NUMERIC(1)); nprotect++;
    PROTECT(xvec = NEW_NUMERIC(nvars)); nprotect++;
    PROTECT(pvec = NEW_NUMERIC(npars)); nprotect++;
    SET_NAMES(xvec,Snames);
    SET_NAMES(pvec,Pnames);

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

    // get function's environment
    PROTECT(rho = (CLOENV(fn))); nprotect++;

    break;

  case 1:			// native code

    // construct state, parameter, covariate indices
    sidx = INTEGER(PROTECT(matchnames(Snames,statenames))); nprotect++;
    pidx = INTEGER(PROTECT(matchnames(Pnames,paramnames))); nprotect++;
    cidx = INTEGER(PROTECT(matchnames(Cnames,covarnames))); nprotect++;

    ff = (pomp_onestep_sim *) R_ExternalPtrAddr(fn);

    break;

  default:
    error("unrecognized 'mode' in 'euler_simulator'");
    break;
  }

  // create array to hold results
  {
    int dim[3] = {nvars, nreps, ntimes};
    PROTECT(X = makearray(3,dim)); nprotect++;
    setrownames(X,Snames,3);
  }

  // copy the start values into the result array
  memcpy(REAL(X),REAL(xstart),nvars*nreps*sizeof(double));

  if (mode==1) {
    set_pomp_userdata(args);
    GetRNGstate();
  }

  // now do computations
  {
    int first = 1;
    int use_names = 0;
    int *posn;
    double *time = REAL(times);
    double *xs = REAL(X);
    double *xt = REAL(X)+nvars*nreps;
    double *cp = REAL(cvec);
    double *ps = REAL(params);
    double t = time[0];
    double dt;
    double *pm, *xm;
    int i, j, k, step;

    for (step = 1; step < ntimes; step++, xs = xt, xt += nvars*nreps) {

      R_CheckUserInterrupt();
	
      if (t > time[step]) {
	error("'times' is not an increasing sequence");
      }

      memcpy(xt,xs,nreps*nvars*sizeof(double));
	
      // set accumulator variables to zero 
      for (j = 0; j < nreps; j++)
	for (i = 0; i < nzeros; i++) 
	  xt[zidx[i]+nvars*j] = 0.0;

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
	error("unrecognized 'method' in 'euler_model_simulator'");
	break;
      }

      for (k = 0; k < nstep; k++) { // loop over Euler steps

	// interpolate the covar functions for the covariates
	table_lookup(&covariate_table,t,cp,0);

	for (j = 0, pm = ps, xm = xt; j < nreps; j++, pm += npars, xm += nvars) { // loop over replicates
	  
	  switch (mode) {

	  case 0: 		// R function

	    {
	      double *xp = REAL(xvec);
	      double *pp = REAL(pvec);
	      double *tp = REAL(tvec);
	      double *dtp = REAL(dtvec);
	      double *ap;
	      
	      *tp = t;
	      *dtp = dt;
	      memcpy(xp,xm,nvars*sizeof(double));
	      memcpy(pp,pm,npars*sizeof(double));
	      
	      if (first) {

	      	PROTECT(ans = eval(fcall,rho));	nprotect++; // evaluate the call
	      	if (LENGTH(ans) != nvars) {
	      	  error("user 'step.fun' returns a vector of %d states but %d are expected: compare initial conditions?",
	      		LENGTH(ans),nvars);
	      	}
		
	      	PROTECT(nm = GET_NAMES(ans)); nprotect++;
	      	use_names = !isNull(nm);
	      	if (use_names) {
	      	  posn = INTEGER(PROTECT(matchnames(Snames,nm))); nprotect++;
	      	}

	      	ap = REAL(AS_NUMERIC(ans));
		
	      	first = 0;

	      } else {
	      
		ap = REAL(AS_NUMERIC(eval(fcall,rho)));

	      }
	      
	      if (use_names) {
	      	for (i = 0; i < nvars; i++) xm[posn[i]] = ap[i];
	      } else {
	      	for (i = 0; i < nvars; i++) xm[i] = ap[i];
	      }

	    }

	    break;
	      
	  case 1: 		// native code

	    (*ff)(xm,pm,sidx,pidx,cidx,ncovars,cp,t,dt);
	    break;

	  default:
	    error("unrecognized 'mode' in 'euler_simulator'");
	    break;
	  }

	}

	t += dt;
	
	if ((method == 0) && (k == nstep-2)) { // penultimate step
	  dt = time[step]-t;
	  t = time[step]-dt;
	}
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
SEXP euler_model_density (SEXP func, 
			  SEXP x, SEXP times, SEXP params, 
			  SEXP statenames, SEXP paramnames, SEXP covarnames,
			  SEXP tcovar, SEXP covar, SEXP log, SEXP args, SEXP gnsi) 
{
  int nprotect = 0;
  int mode;
  int give_log;
  int nvars, npars, nreps, ntimes, ncovars, covlen;
  pomp_onestep_pdf *ff = NULL;
  SEXP t1vec, t2vec, x1vec, x2vec, pvec, cvec;
  SEXP Snames, Pnames, Cnames;
  SEXP rho, fcall, fn;
  SEXP F;
  int *pidx = 0, *sidx = 0, *cidx = 0;

  {
    int *dim;
    dim = INTEGER(GET_DIM(x)); nvars = dim[0]; nreps = dim[1];
    dim = INTEGER(GET_DIM(params)); npars = dim[0];
    dim = INTEGER(GET_DIM(covar)); covlen = dim[0]; ncovars = dim[1];
    ntimes = LENGTH(times);
  }

  PROTECT(Snames = GET_ROWNAMES(GET_DIMNAMES(x))); nprotect++;
  PROTECT(Pnames = GET_ROWNAMES(GET_DIMNAMES(params))); nprotect++;
  PROTECT(Cnames = GET_COLNAMES(GET_DIMNAMES(covar))); nprotect++;

  // set up the covariate table
  struct lookup_table covariate_table = {covlen, ncovars, 0, REAL(tcovar), REAL(covar)};

  // vector for interpolated covariates
  PROTECT(cvec = NEW_NUMERIC(ncovars)); nprotect++;
  SET_NAMES(cvec,Cnames);

  PROTECT(fn = pomp_fun_handler(func,gnsi,&mode)); nprotect++;

  give_log = *(INTEGER(log));

  switch (mode) {

  case 0:			// R function

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

  case 1:			// native code

    // construct state, parameter, covariate indices
    sidx = INTEGER(PROTECT(matchnames(Snames,statenames))); nprotect++;
    pidx = INTEGER(PROTECT(matchnames(Pnames,paramnames))); nprotect++;
    cidx = INTEGER(PROTECT(matchnames(Cnames,covarnames))); nprotect++;

    ff = (pomp_onestep_pdf *) R_ExternalPtrAddr(fn);

    break;

  default:
    error("unrecognized 'mode' in 'euler_model_density'");
    break;
  }

  // create array to hold results
  {
    int dim[2] = {nreps, ntimes-1};
    PROTECT(F = makearray(2,dim)); nprotect++;
  }

  switch (mode) {

  case 0:			// R function

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
	table_lookup(&covariate_table,*t1p,cp,0);
    
	for (j = 0, ps = REAL(params); j < nreps; j++, fs++, x1s += nvars, x2s += nvars, ps += npars) { // loop over replicates
      
	  memcpy(x1p,x1s,nvars*sizeof(double));
	  memcpy(x2p,x2s,nvars*sizeof(double));
	  memcpy(pp,ps,npars*sizeof(double));

	  *fs = *(REAL(AS_NUMERIC(eval(fcall,rho))));
      
	  if (!give_log) *fs = exp(*fs);
      
	}
      }
    }

    break;

  case 1:			// native code

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
	table_lookup(&covariate_table,*t1s,cp,0);
    
	for (j = 0, ps = REAL(params); j < nreps; j++, fs++, x1s += nvars, x2s += nvars, ps += npars) { // loop over replicates
      
	  (*ff)(fs,x1s,x2s,*t1s,*t2s,ps,sidx,pidx,cidx,ncovars,cp);

	  if (!give_log) *fs = exp(*fs);
      
	}
      }
    }

    unset_pomp_userdata();

    break;

  default:
    error("unrecognized 'mode' in 'euler_model_density'");
    break;

  }

  UNPROTECT(nprotect);
  return F;
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

