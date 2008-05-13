// dear emacs, please treat this as -*- C++ -*-

#include "pomp.h"

// take nstep Euler-Poisson steps from t1 to t2
static void euler_simulator (euler_step_sim *estep,
			     double *x, double *xstart, double *times, double *params, 
			     int *ndim, double *deltat,
			     int *stateindex, int *parindex, int *zeroindex,
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
  double *covar[covdim];
  int k, p, step;
  double dt;

  GetRNGstate();		// initialize R's pseudorandom number generator

  for (k = 0; k < covdim; k++) covar[k] = &covar_table[k*covlen]; // set up the covariate table
  struct lookup_table covariate_table = {covlen, covdim, 0, time_table, covar};
  
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

	(*estep)(xp,pp,stateindex,parindex,covdim,covar_fn,t,dt);

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

SEXP euler_model_simulator (SEXP func, 
			    SEXP xstart, SEXP times, SEXP params, 
			    SEXP dt, 
			    SEXP statenames, SEXP paramnames, SEXP zeronames,
			    SEXP tbasis, SEXP basis) 
{
  int nprotect = 0;
  int *dim, xdim[3], ndim[7];
  int nvar, npar, nrep, ntimes;
  int poplen, baslen, basdim;
  int nstates = length(statenames);
  int nparams = length(paramnames);
  int nzeros = length(zeronames);
  euler_step_sim *ff;
  SEXP X, dimX, pnm, snm, pindex, sindex, zindex;
  int k;

  if (inherits(func,"NativeSymbol")) {
    ff = (euler_step_sim *) R_ExternalPtrAddr(func);
  } else {
    error("illegal input: supplied function is not a compiled function");
  }

  dim = INTEGER(GET_DIM(xstart)); nvar = dim[0]; nrep = dim[1];
  dim = INTEGER(GET_DIM(params)); npar = dim[0];
  dim = INTEGER(GET_DIM(basis)); baslen = dim[0]; basdim = dim[1];
  ntimes = length(times);
  xdim[0] = nvar; xdim[1] = nrep; xdim[2] = ntimes;
  PROTECT(X = makearray(3,xdim)); nprotect++;
  setrownames(X,GET_ROWNAMES(GET_DIMNAMES(xstart)),3);
  PROTECT(sindex = matchrownames(xstart,statenames)); nprotect++;
  PROTECT(pindex = matchrownames(params,paramnames)); nprotect++;
  PROTECT(zindex = matchrownames(xstart,zeronames)); nprotect++;
  ndim[0] = nvar; ndim[1] = npar; ndim[2] = nrep; ndim[3] = ntimes; 
  ndim[4] = baslen; ndim[5] = basdim; ndim[6] = nzeros;
  euler_simulator(ff,REAL(X),REAL(xstart),REAL(times),REAL(params),
		  ndim,REAL(dt),
		  INTEGER(sindex),INTEGER(pindex),INTEGER(zindex),
		  REAL(tbasis),REAL(basis));
  UNPROTECT(nprotect);
  return X;
}

// take nstep Euler-Poisson steps from t1 to t2
static void euler_densities (euler_step_pdf *estep,
			     double *f, 
			     double *x, double *times, double *params, 
			     int *ndim,
			     int *stateindex, int *parindex,
			     double *time_table, double *covar_table,
			     int *give_log)
{
  double t, *x1p, *x2p, *pp, *fp;
  int nvar = ndim[0];
  int npar = ndim[1];
  int nrep = ndim[2];
  int ntimes = ndim[3];
  int covlen = ndim[4];
  int covdim = ndim[5];
  double covar_fn[covdim];
  double *covar[covdim];
  int k, p, step;
  double dt;

  // set up the covariate table
  for (k = 0; k < covdim; k++) covar[k] = &covar_table[k*covlen]; 
  struct lookup_table covariate_table = {covlen, covdim, 0, time_table, covar};
  
  for (step = 0; step < ntimes-1; step++) { // loop over times

    R_CheckUserInterrupt();	// check for user interrupt

    t = times[step];
    dt = times[step+1]-t;

    // interpolate the covar functions for the covariates
    if (covdim > 0) 
      table_lookup(&covariate_table,t,covar_fn,0);
    
    for (p = 0; p < nrep; p++) { // loop over replicates
      
      fp = &f[p+nrep*step];
      x1p = &x[nvar*(p+nrep*step)];
      x2p = &x[nvar*(p+nrep*(step+1))];
      pp = &params[npar*p];
      
      (*estep)(fp,x1p,x2p,pp,
	       stateindex,parindex,
	       covdim,covar_fn,t,dt);
      
      if (!(*give_log)) *fp = exp(*fp);
      
    }
  }
}

SEXP euler_model_density (SEXP func, 
			  SEXP x, SEXP times, SEXP params, 
			  SEXP statenames, SEXP paramnames,
			  SEXP tbasis, SEXP basis, SEXP log) 
{
  int nprotect = 0;
  int *dim, fdim[2], ndim[6];
  int nvar, npar, nrep, ntimes;
  int baslen, basdim;
  int nstates = length(statenames);
  int nparams = length(paramnames);
  euler_step_pdf *ff;
  SEXP F, dimF, pnm, snm, pindex, sindex;
  int k;

  if (inherits(func,"NativeSymbol")) {
    ff = (euler_step_pdf *) R_ExternalPtrAddr(func);
  } else {
    error("illegal input: supplied function is not a compiled function");
  }

  dim = INTEGER(GET_DIM(x)); nvar = dim[0]; nrep = dim[1];
  dim = INTEGER(GET_DIM(params)); npar = dim[0];
  dim = INTEGER(GET_DIM(basis)); baslen = dim[0]; basdim = dim[1];
  ntimes = length(times);
  fdim[0] = nrep; fdim[1] = ntimes-1;
  PROTECT(F = makearray(2,fdim)); nprotect++;
  PROTECT(sindex = matchrownames(x,statenames)); nprotect++;
  PROTECT(pindex = matchrownames(params,paramnames)); nprotect++;
  ndim[0] = nvar; ndim[1] = npar; ndim[2] = nrep; ndim[3] = ntimes; 
  ndim[4] = baslen; ndim[5] = basdim;
  euler_densities(ff,REAL(F),REAL(x),REAL(times),REAL(params),
		  ndim,INTEGER(sindex),INTEGER(pindex),
		  REAL(tbasis),REAL(basis),INTEGER(log));
  UNPROTECT(nprotect);
  return F;
}
