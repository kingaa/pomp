// dear emacs, please treat this as -*- C++ -*-

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

#include "pomp_internal.h"

static void dens_meas (pomp_measure_model_density *f,
		       double *lik, double *y, 
		       double *x, double *times, double *params, 
		       int *give_log, int *ndim,
		       int *obsindex, int *stateindex, int *parindex, int *covindex,
		       double *time_table, double *covar_table)
{
  double t, *lp, *xp, *pp, *yp;
  int nvar = ndim[0];
  int npar = ndim[1];
  int nrepp = ndim[2];
  int nrepx = ndim[3];
  int ntimes = ndim[4];
  int covlen = ndim[5];
  int covdim = ndim[6];
  int nobs = ndim[7];
  double covar_fn[covdim];
  int nrep;
  int k, p;
  
  // set up the covariate table
  struct lookup_table covariate_table = {covlen, covdim, 0, time_table, covar_table};

  nrep = (nrepp > nrepx) ? nrepp : nrepx;
  
  for (k = 0; k < ntimes; k++) { // loop over times

    R_CheckUserInterrupt();	// check for user interrupt

    t = times[k];
    yp = &y[nobs*k];

    // interpolate the covar functions for the covariates
    if (covdim > 0) 
      table_lookup(&covariate_table,t,covar_fn,0);
    
    for (p = 0; p < nrep; p++) { // loop over replicates
      
      lp = &lik[p+nrep*k];
      xp = &x[nvar*((p%nrepx)+nrepx*k)];
      pp = &params[npar*(p%nrepp)];

      (*f)(lp,yp,xp,pp,*give_log,obsindex,stateindex,parindex,covindex,covdim,covar_fn,t);
      
    }
  }
}

// these global objects will pass the needed information to the user-defined function
// each of these is allocated once, globally, and refilled many times
static SEXP _pomp_dmeas_Yvec;	// observable vector
static SEXP _pomp_dmeas_Xvec;	// state variable vector
static SEXP _pomp_dmeas_Pvec;	// parameter vector
static SEXP _pomp_dmeas_Cvec;	// covariate vector
static SEXP _pomp_dmeas_time;	// time
static int  _pomp_dmeas_nvar;	// number of state variables
static int  _pomp_dmeas_npar;	// number of parameters
static int  _pomp_dmeas_nobs;	// number of observables
static SEXP _pomp_dmeas_envir;	// function's environment
static SEXP _pomp_dmeas_fcall;	// function call

#define YVEC    (_pomp_dmeas_Yvec)
#define XVEC    (_pomp_dmeas_Xvec)
#define PVEC    (_pomp_dmeas_Pvec)
#define CVEC    (_pomp_dmeas_Cvec)
#define TIME    (_pomp_dmeas_time)
#define NVAR    (_pomp_dmeas_nvar)
#define NPAR    (_pomp_dmeas_npar)
#define NOBS    (_pomp_dmeas_nobs)
#define RHO     (_pomp_dmeas_envir)
#define FCALL   (_pomp_dmeas_fcall)

// this is the measurement pdf evaluated when the user supplies an R function
// (and not a native routine)
// Note that obsindex, stateindex, parindex, covindex are ignored.
static void default_meas_dens (double *lik, double *y, double *x, double *p, int give_log,
			       int *obsindex, int *stateindex, int *parindex, int *covindex,
			       int covdim, double *covar, double t)
{
  int nprotect = 0;
  int k;
  double *xp;
  SEXP ans;
  xp = REAL(YVEC);
  for (k = 0; k < NOBS; k++) xp[k] = y[k];
  xp = REAL(XVEC);
  for (k = 0; k < NVAR; k++) xp[k] = x[k];
  xp = REAL(PVEC);
  for (k = 0; k < NPAR; k++) xp[k] = p[k];
  xp = REAL(CVEC);
  for (k = 0; k < covdim; k++) xp[k] = covar[k];
  xp = REAL(TIME);
  xp[0] = t;
  PROTECT(ans = AS_NUMERIC(eval(FCALL,RHO))); nprotect++; // evaluate the call
  if (LENGTH(ans) != 1)
    error("dmeasure error: user 'dmeasure' must return a scalar");
  xp = REAL(ans);
  *lik = *xp;
  UNPROTECT(nprotect);
}

SEXP do_dmeasure (SEXP object, SEXP y, SEXP x, SEXP times, SEXP params, SEXP log, SEXP fun)
{
  int nprotect = 0;
  int *dim, nvars, npars, nrepsp, nrepsx, nreps, ntimes, covlen, covdim, nobs;
  int ndim[7];
  SEXP F, fn;
  SEXP tcovar, covar;
  SEXP statenames, paramnames, covarnames, obsnames;
  SEXP sindex, pindex, cindex, oindex;
  int *sidx, *pidx, *cidx, *oidx;
  SEXP Xnames, Ynames, Pnames, Cnames;
  int use_native;
  int nstates, nparams, ncovars, nobsers;
  pomp_measure_model_density *ff = NULL;

  PROTECT(times = AS_NUMERIC(times)); nprotect++;
  ntimes = LENGTH(times);
  if (ntimes < 1)
    error("dmeasure error: no work to do");

  PROTECT(y = as_matrix(y)); nprotect++;
  dim = INTEGER(GET_DIM(y));
  nobs = dim[0];

  if (ntimes != dim[1])
    error("dmeasure error: length of 'times' and 2nd dimension of 'y' do not agree");

  PROTECT(x = as_state_array(x)); nprotect++;
  dim = INTEGER(GET_DIM(x));
  nvars = dim[0]; nrepsx = dim[1]; 

  if (ntimes != dim[2])
    error("dmeasure error: length of 'times' and 3rd dimension of 'x' do not agree");

  PROTECT(params = as_matrix(params)); nprotect++;
  dim = INTEGER(GET_DIM(params));
  npars = dim[0]; nrepsp = dim[1]; 

  nreps = (nrepsp > nrepsx) ? nrepsp : nrepsx;

  if ((nreps % nrepsp != 0) || (nreps % nrepsx != 0))
    error("dmeasure error: larger number of replicates is not a multiple of smaller");

  PROTECT(tcovar =  GET_SLOT(object,install("tcovar"))); nprotect++;
  PROTECT(covar =  GET_SLOT(object,install("covar"))); nprotect++;
  PROTECT(obsnames =  GET_SLOT(object,install("obsnames"))); nprotect++;
  PROTECT(statenames =  GET_SLOT(object,install("statenames"))); nprotect++;
  PROTECT(paramnames =  GET_SLOT(object,install("paramnames"))); nprotect++;
  PROTECT(covarnames =  GET_SLOT(object,install("covarnames"))); nprotect++;
  nobsers = LENGTH(obsnames);
  nstates = LENGTH(statenames);
  nparams = LENGTH(paramnames);
  ncovars = LENGTH(covarnames);

  dim = INTEGER(GET_DIM(covar)); covlen = dim[0]; covdim = dim[1];
  PROTECT(Ynames = GET_ROWNAMES(GET_DIMNAMES(y))); nprotect++;
  PROTECT(Xnames = GET_ROWNAMES(GET_DIMNAMES(x))); nprotect++;
  PROTECT(Pnames = GET_ROWNAMES(GET_DIMNAMES(params))); nprotect++;
  PROTECT(Cnames = GET_COLNAMES(GET_DIMNAMES(covar))); nprotect++;

  PROTECT(fn = VECTOR_ELT(fun,0)); nprotect++;
  use_native = INTEGER(VECTOR_ELT(fun,1))[0];

  switch (use_native) {
  case 0:			// use R function
    ff = (pomp_measure_model_density *) default_meas_dens;
    PROTECT(RHO = (CLOENV(fn))); nprotect++;
    NVAR = nvars;			// for internal use
    NPAR = npars;			// for internal use
    NOBS = nobs;			// for internal use
    PROTECT(TIME = NEW_NUMERIC(1)); nprotect++;	// for internal use
    PROTECT(YVEC = NEW_NUMERIC(nobs)); nprotect++; // for internal use
    PROTECT(XVEC = NEW_NUMERIC(nvars)); nprotect++; // for internal use
    PROTECT(PVEC = NEW_NUMERIC(npars)); nprotect++; // for internal use
    PROTECT(CVEC = NEW_NUMERIC(covdim)); nprotect++; // for internal use
    SET_NAMES(YVEC,Ynames); // make sure the names attribute is copied
    SET_NAMES(XVEC,Xnames); // make sure the names attribute is copied
    SET_NAMES(PVEC,Pnames); // make sure the names attribute is copied
    SET_NAMES(CVEC,Cnames); // make sure the names attribute is copied
    // set up the function call
    PROTECT(FCALL = VectorToPairList(GET_SLOT(object,install("userdata")))); nprotect++;
    PROTECT(FCALL = LCONS(CVEC,FCALL)); nprotect++;
    SET_TAG(FCALL,install("covars"));
    PROTECT(FCALL = LCONS(AS_LOGICAL(log),FCALL)); nprotect++;
    SET_TAG(FCALL,install("log"));
    PROTECT(FCALL = LCONS(PVEC,FCALL)); nprotect++;
    SET_TAG(FCALL,install("params"));
    PROTECT(FCALL = LCONS(TIME,FCALL)); nprotect++;
    SET_TAG(FCALL,install("t"));
    PROTECT(FCALL = LCONS(XVEC,FCALL)); nprotect++;
    SET_TAG(FCALL,install("x"));
    PROTECT(FCALL = LCONS(YVEC,FCALL)); nprotect++;
    SET_TAG(FCALL,install("y"));
    PROTECT(FCALL = LCONS(fn,FCALL)); nprotect++;
    break;
  case 1:			// use native routine
    ff = (pomp_measure_model_density *) R_ExternalPtrAddr(fn);
    break;
  default:
    error("unrecognized 'use' slot in 'dmeasure'");
    break;
  }

  ndim[0] = nreps; ndim[1] = ntimes;
  PROTECT(F = makearray(2,ndim)); nprotect++; 

  if (nobsers > 0) {
    PROTECT(oindex = matchnames(Ynames,obsnames)); nprotect++;
    oidx = INTEGER(oindex);
  } else {
    oidx = 0;
  }

  if (nstates > 0) {
    PROTECT(sindex = matchnames(Xnames,statenames)); nprotect++;
    sidx = INTEGER(sindex);
  } else {
    sidx = 0;
  }

  if (nparams > 0) {
    PROTECT(pindex = matchnames(Pnames,paramnames)); nprotect++;
    pidx = INTEGER(pindex);
  } else {
    pidx = 0;
  }
  
  if (ncovars > 0) {
    PROTECT(cindex = matchnames(Cnames,covarnames)); nprotect++;
    cidx = INTEGER(cindex);
  } else {
    cidx = 0;
  }
  
  ndim[0] = nvars; ndim[1] = npars; ndim[2] = nrepsp; ndim[3] = nrepsx; ndim[4] = ntimes; 
  ndim[5] = covlen; ndim[6] = covdim; ndim[7] = nobs;

  if (use_native) {
    PROTECT(FCALL = VectorToPairList(GET_SLOT(object,install("userdata")))); nprotect++;
    set_pomp_userdata(FCALL);
  }

  dens_meas(ff,REAL(F),REAL(y),REAL(x),REAL(times),REAL(params),INTEGER(log),ndim,
	    oidx,sidx,pidx,cidx,REAL(tcovar),REAL(covar));
  
  if (use_native) {
    unset_pomp_userdata();
  }

  UNPROTECT(nprotect);
  return F;
}

#undef YVEC
#undef XVEC
#undef PVEC
#undef CVEC
#undef TIME
#undef NVAR
#undef NPAR
#undef NOBS
#undef RHO
#undef FCALL
