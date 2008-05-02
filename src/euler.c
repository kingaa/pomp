// dear emacs, please treat this as -*- C++ -*-

#include "euler.h"
#include "interp.h"

static SEXP makearray(int rank, int *dim);
static SEXP matchrownames (SEXP x, SEXP names);
static void setrownames (SEXP x, SEXP names);
static void euler_simulator (euler_step *estep,
			     double *x, double *xstart, double *times, double *params, 
			    int *ndim, double *deltat,
			    int *stateindex, int *parindex, int *zeroindex,
			    double *time_table, double *covar_table);

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
  euler_step *ff;
  SEXP X, dimX, pnm, snm, pindex, sindex, zindex;
  int k;

  if (inherits(func,"NativeSymbol")) {
    ff = (euler_step *) R_ExternalPtrAddr(func);
  } else {
    error("illegal input: supplied function is not a compiled function");
  }

  dim = INTEGER(GET_DIM(xstart)); nvar = dim[0]; nrep = dim[1];
  dim = INTEGER(GET_DIM(params)); npar = dim[0];
  dim = INTEGER(GET_DIM(basis)); baslen = dim[0]; basdim = dim[1];
  ntimes = length(times);
  xdim[0] = nvar; xdim[1] = nrep; xdim[2] = ntimes;
  PROTECT(X = makearray(3,xdim)); nprotect++;
  setrownames(X,GET_ROWNAMES(GET_DIMNAMES(xstart)));
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

static SEXP makearray (int rank, int *dim) {
  int nprotect = 0;
  int *dimp, k;
  SEXP dimx, x;
  PROTECT(dimx = NEW_INTEGER(rank)); nprotect++;
  dimp = INTEGER(dimx); 
  for (k = 0; k < rank; k++) dimp[k] = dim[k];
  PROTECT(x = allocArray(REALSXP,dimx)); nprotect++;
  UNPROTECT(nprotect);
  return x;
}

static SEXP matchrownames (SEXP x, SEXP names) {
  int nprotect = 0;
  int n = length(names);
  int *idx, k;
  SEXP index, nm;
  PROTECT(nm = AS_CHARACTER(names)); nprotect++;
  PROTECT(index = match(GET_ROWNAMES(GET_DIMNAMES(x)),names,0)); nprotect++;
  idx = INTEGER(index);
  for (k = 0; k < n; k++) {
    if (idx[k]==0) error("variable %s not specified",STRING_ELT(nm,k));
    idx[k] -= 1;
  }
  UNPROTECT(nprotect);
  return index;
}

static void setrownames (SEXP x, SEXP names) {
  int nprotect = 0;
  int n = length(names);
  int k;
  SEXP dimnms, nm;
  PROTECT(nm = AS_CHARACTER(names)); nprotect++;
  PROTECT(dimnms = allocVector(VECSXP,3)); nprotect++;
  SET_ELEMENT(dimnms,0,nm);	// set row names
  SET_DIMNAMES(x,dimnms);
  UNPROTECT(nprotect);
}

// take nstep Euler-Poisson steps from t1 to t2
static void euler_simulator (euler_step *estep,
			     double *x, double *xstart, double *times, double *params, 
			     int *ndim, double *deltat,
			     int *stateindex, int *parindex, int *zeroindex,
			     double *time_table, double *covar_table)
{
  double pop, dpopdt, beta;
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

