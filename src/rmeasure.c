// dear emacs, please treat this as -*- C++ -*-

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

#include "pomp_internal.h"

static void simul_meas (pomp_measure_model_simulator *f,
			double *y, 
			double *x, double *times, double *params, 
			int *ndim,
			int *stateindex, int *parindex, int *covindex,
			double *time_table, double *covar_table)
{
  double t, *xp, *pp, *yp;
  int nvar = ndim[0];
  int npar = ndim[1];
  int nrep = ndim[2];
  int ntimes = ndim[3];
  int covlen = ndim[4];
  int covdim = ndim[5];
  int nobs = ndim[6];
  double covar_fn[covdim];
  int k, p;
  
  // set up the covariate table
  struct lookup_table covariate_table = {covlen, covdim, 0, time_table, covar_table};
  
  for (k = 0; k < ntimes; k++) { // loop over times

    R_CheckUserInterrupt();	// check for user interrupt

    t = times[k];

    // interpolate the covar functions for the covariates
    if (covdim > 0) 
      table_lookup(&covariate_table,t,covar_fn,0);
    
    for (p = 0; p < nrep; p++) { // loop over replicates
      
      yp = &y[nobs*(p+nrep*k)];
      xp = &x[nvar*(p+nrep*k)];
      pp = &params[npar*p];

      (*f)(yp,xp,pp,stateindex,parindex,covindex,covdim,covar_fn,t);
      
    }
  }
}

// these global objects will pass the needed information to the user-defined function
// each of these is allocated once, globally, and refilled many times
static SEXP _pomp_rmeas_Xvec;	// state variable vector
static SEXP _pomp_rmeas_Pvec;	// parameter vector
static SEXP _pomp_rmeas_Cvec;	// covariate vector
static SEXP _pomp_rmeas_time;	// time
static int  _pomp_rmeas_nvar;	// number of state variables
static int  _pomp_rmeas_npar;	// number of parameters
static int  _pomp_rmeas_nobs;	// number of observables
static int  _pomp_rmeas_first;	// first evaluation?
static SEXP _pomp_rmeas_onames;	// names of observable
static int *_pomp_rmeas_oindex;	// indices of observables
static SEXP _pomp_rmeas_envir;	// function's environment
static SEXP _pomp_rmeas_fcall;	// function call

#define XVEC    (_pomp_rmeas_Xvec)
#define PVEC    (_pomp_rmeas_Pvec)
#define CVEC    (_pomp_rmeas_Cvec)
#define TIME    (_pomp_rmeas_time)
#define NVAR    (_pomp_rmeas_nvar)
#define NPAR    (_pomp_rmeas_npar)
#define NOBS    (_pomp_rmeas_nobs)
#define FIRST   (_pomp_rmeas_first)
#define OBNM    (_pomp_rmeas_onames)
#define OIDX    (_pomp_rmeas_oindex)
#define RHO     (_pomp_rmeas_envir)
#define FCALL   (_pomp_rmeas_fcall)

// this is the measurement simulator evaluated when the user supplies an R function
// (and not a native routine)
// Note that stateindex, parindex, covindex are ignored.
static void default_meas_sim (double *y, double *x, double *p, 
			      int *stateindex, int *parindex, int *covindex,
			      int covdim, double *covar, double t)
{
  int nprotect = 0;
  int k, *op;
  double *xp;
  SEXP ans, nm, oidx;
  xp = REAL(XVEC);
  for (k = 0; k < NVAR; k++) xp[k] = x[k];
  xp = REAL(PVEC);
  for (k = 0; k < NPAR; k++) xp[k] = p[k];
  xp = REAL(CVEC);
  for (k = 0; k < covdim; k++) xp[k] = covar[k];
  xp = REAL(TIME);
  xp[0] = t;

  PROTECT(ans = eval(FCALL,RHO)); nprotect++; // evaluate the call

  if (FIRST) {
    if (LENGTH(ans) != NOBS)
      error("rmeasure error: user 'rmeasure' must return a vector of length ",NOBS);
    PROTECT(nm = GET_NAMES(ans)); nprotect++;
    if (!isNull(nm)) {		// match names against names from data slot
      PROTECT(oidx = matchnames(OBNM,nm)); nprotect++;
      op = INTEGER(oidx);
      for (k = 0; k < NOBS; k++) OIDX[k] = op[k];
    } else {
      for (k = 0; k < NOBS; k++) OIDX[k] = k;
    }
    FIRST = 0;
  }

  xp = REAL(AS_NUMERIC(ans));
  for (k = 0; k < NOBS; k++) y[OIDX[k]] = xp[k];
  UNPROTECT(nprotect);
}

SEXP do_rmeasure (SEXP object, SEXP x, SEXP times, SEXP params)
{
  int nprotect = 0;
  int *dim, nvars, npars, nreps, ntimes, covlen, covdim, nobs;
  int ndim[7];
  SEXP Y, fn;
  SEXP dimP, dimX, dimD;
  SEXP tcovar, covar;
  SEXP statenames, paramnames, covarnames;
  SEXP sindex, pindex, cindex;
  int *sidx, *pidx, *cidx;
  SEXP Xnames, Pnames, Cnames;
  int use_native;
  int nstates, nparams, ncovars;
  pomp_measure_model_simulator *ff = NULL;

  PROTECT(times = AS_NUMERIC(times)); nprotect++;
  ntimes = length(times);
  if (ntimes < 1)
    error("rmeasure error: no work to do");

  PROTECT(dimX = GET_DIM(x)); nprotect++;
  if ((isNull(dimX)) || (length(dimX)!=3))
    error("rmeasure error: 'x' must be a rank-3 array");
  dim = INTEGER(dimX); nvars = dim[0]; nreps = dim[1];
  if (ntimes != dim[2])
    error("rprocess error: length of 'times' and 3rd dimension of 'x' do not agree");

  PROTECT(dimP = GET_DIM(params)); nprotect++;
  if ((isNull(dimP)) || (length(dimP)!=2))
    error("rmeasure error: 'params' must be a rank-2 array");
  dim = INTEGER(dimP); npars = dim[0];
  if (nreps != dim[1])
    error("rmeasure error: 2nd dimensions of 'params' and 'x' do not agree");

  PROTECT(dimD = GET_DIM(GET_SLOT(object,install("data")))); nprotect++;
  dim = INTEGER(dimD); nobs = dim[0];

  PROTECT(tcovar =  GET_SLOT(object,install("tcovar"))); nprotect++;
  PROTECT(covar =  GET_SLOT(object,install("covar"))); nprotect++;
  PROTECT(Cnames = GET_COLNAMES(GET_DIMNAMES(covar))); nprotect++;
  PROTECT(statenames =  GET_SLOT(object,install("statenames"))); nprotect++;
  PROTECT(paramnames =  GET_SLOT(object,install("paramnames"))); nprotect++;
  PROTECT(covarnames =  GET_SLOT(object,install("covarnames"))); nprotect++;
  nstates = LENGTH(statenames);
  nparams = LENGTH(paramnames);
  ncovars = LENGTH(covarnames);

  dim = INTEGER(GET_DIM(covar)); covlen = dim[0]; covdim = dim[1];
  PROTECT(Xnames = GET_ROWNAMES(GET_DIMNAMES(x))); nprotect++;
  PROTECT(Pnames = GET_ROWNAMES(GET_DIMNAMES(params))); nprotect++;
  PROTECT(Cnames = GET_COLNAMES(GET_DIMNAMES(covar))); nprotect++;
  PROTECT(OBNM = GET_ROWNAMES(GET_DIMNAMES(GET_SLOT(object,install("data"))))); nprotect++;

  PROTECT(
	  fn = pomp_fun_handler(
				GET_SLOT(object,install("rmeasure")),
				&use_native
				)
	  ); nprotect++;

  if (use_native) {		// use native routine
    ff = (pomp_measure_model_simulator *) R_ExternalPtrAddr(fn);
    OIDX = 0;
  } else {			// use R function
    ff = (pomp_measure_model_simulator *) default_meas_sim;
    PROTECT(RHO = (CLOENV(fn))); nprotect++;
    NVAR = nvars;			// for internal use
    NPAR = npars;			// for internal use
    NOBS = nobs;			// for internal use
    PROTECT(TIME = NEW_NUMERIC(1)); nprotect++;	// for internal use
    PROTECT(XVEC = NEW_NUMERIC(nvars)); nprotect++; // for internal use
    PROTECT(PVEC = NEW_NUMERIC(npars)); nprotect++; // for internal use
    PROTECT(CVEC = NEW_NUMERIC(covdim)); nprotect++; // for internal use
    SET_NAMES(XVEC,Xnames); // make sure the names attribute is copied
    SET_NAMES(PVEC,Pnames); // make sure the names attribute is copied
    SET_NAMES(CVEC,Cnames); // make sure the names attribute is copied
    // set up the function call
    PROTECT(FCALL = VectorToPairList(GET_SLOT(object,install("userdata")))); nprotect++;
    PROTECT(FCALL = LCONS(CVEC,FCALL)); nprotect++;
    SET_TAG(FCALL,install("covars"));
    PROTECT(FCALL = LCONS(PVEC,FCALL)); nprotect++;
    SET_TAG(FCALL,install("params"));
    PROTECT(FCALL = LCONS(TIME,FCALL)); nprotect++;
    SET_TAG(FCALL,install("t"));
    PROTECT(FCALL = LCONS(XVEC,FCALL)); nprotect++;
    SET_TAG(FCALL,install("x"));
    PROTECT(FCALL = LCONS(fn,FCALL)); nprotect++;
    OIDX = (int *) Calloc(nobs,int);
    FIRST = 1;
  }

  ndim[0] = nobs; ndim[1] = nreps; ndim[2] = ntimes;
  PROTECT(Y = makearray(3,ndim)); nprotect++; 
  setrownames(Y,OBNM,3);

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

  ndim[0] = nvars; ndim[1] = npars; ndim[2] = nreps; ndim[3] = ntimes; 
  ndim[4] = covlen; ndim[5] = covdim; ndim[6] = nobs;

  simul_meas(ff,REAL(Y),REAL(x),REAL(times),REAL(params),ndim,
	     sidx,pidx,cidx,
	     REAL(tcovar),REAL(covar));

  if (OIDX != 0) Free(OIDX);

  UNPROTECT(nprotect);
  return Y;
}

#undef FCALL
#undef XVEC
#undef PVEC
#undef CVEC
#undef TIME
#undef NVAR
#undef NPAR
#undef RHO
