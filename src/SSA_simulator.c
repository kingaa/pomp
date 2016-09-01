// dear emacs, please treat this as -*- C++ -*-

#include "pomp_internal.h"
#include <R_ext/Constants.h>

// these global objects will pass the needed information to the user-defined function (see 'default_ssa_internal_fn')
// each of these is allocated once, globally, and refilled many times
static SEXP ssa_internal_jrate;	 // reaction number
static SEXP ssa_internal_Xvec;	 // state variable vector
static SEXP ssa_internal_Pvec;	 // parameter vector
static SEXP ssa_internal_Cvec;	// covariate vector
static SEXP ssa_internal_time;	// time
static int  ssa_internal_nvar;	// number of state variables
static int  ssa_internal_npar;	// number of parameters
static SEXP ssa_internal_envir;	// function's environment
static SEXP ssa_internal_fcall;	// function call
static int  ssa_internal_first;	// first evaluation?
static pomp_ssa_rate_fn *ssa_internal_rxrate; // function computing reaction rates

#define FIRST   (ssa_internal_first)
#define JRATE   (ssa_internal_jrate)
#define XVEC    (ssa_internal_Xvec)
#define PVEC    (ssa_internal_Pvec)
#define CVEC    (ssa_internal_Cvec)
#define TIME    (ssa_internal_time)
#define NVAR    (ssa_internal_nvar)
#define NPAR    (ssa_internal_npar)
#define RHO     (ssa_internal_envir)
#define FCALL   (ssa_internal_fcall)
#define RXR     (ssa_internal_rxrate)

static double default_ssa_rate_fn (int j, double t, const double *x, const double *p,
				   int *stateindex, int *parindex, int *covindex,
				   int ncovar, double *covar)
{
  int nprotect = 0;
  int *op, k;
  double *xp, rate;
  SEXP ans;
  op = INTEGER(JRATE);
  op[0] = j;
  xp = REAL(XVEC);
  for (k = 0; k < NVAR; k++) xp[k] = x[k];
  xp = REAL(TIME);
  xp[0] = t;
  xp = REAL(PVEC);
  for (k = 0; k < NPAR; k++) xp[k] = p[k];
  xp = REAL(CVEC);
  for (k = 0; k < ncovar; k++) xp[k] = covar[k];

  PROTECT(ans = eval(FCALL,RHO)); nprotect++; // evaluate the call

  if (FIRST) {
    if (LENGTH(ans) != 1) {
      UNPROTECT(nprotect);
      errorcall(R_NilValue,"user 'rates' must return a single scalar rate.");
    }
    FIRST = 0;
  }

  rate = *(REAL(AS_NUMERIC(ans)));
  UNPROTECT(nprotect);
  return rate;
}

void SSA (pomp_ssa_rate_fn *ratefun, int irep,
	  int nvar, int nevent, int npar, int nrep, int ntimes,
	  int method,
	  double *xstart, double *times, double *params, double *xout,
	  double *e, double *v, double *d,
	  int ndeps, int *ideps, int nzero, int *izero,
	  int *istate, int *ipar, int ncovar, int *icovar,
	  int lcov, int mcov, double *tcov, double *cov);

SEXP SSA_simulator (SEXP func, SEXP mflag, SEXP xstart, SEXP times, SEXP params,
		    SEXP e, SEXP vmatrix, SEXP dmatrix, SEXP deps, SEXP tcovar, SEXP covar,
		    SEXP zeronames, SEXP args, SEXP gnsi)
{
  int nprotect = 0;
  int *dim, xdim[3];
  int nvar, nevent, npar, nrep, ntimes, ndeps;
  int covlen, covdim;
  SEXP statenames, paramnames, covarnames;
  int nstates, nparams, ncovars;
  int nzeros = LENGTH(zeronames);
  int use_native = 0;
  SEXP X, pindex, sindex, cindex, zindex;
  int *sidx, *pidx, *cidx, *zidx, *didx = 0;
  SEXP fn, Snames, Pnames, Cnames;

  int method = *(INTEGER(mflag));
  dim = INTEGER(GET_DIM(xstart)); nvar = dim[0]; nrep = dim[1];
  dim = INTEGER(GET_DIM(params)); npar = dim[0];
  dim = INTEGER(GET_DIM(covar)); covlen = dim[0]; covdim = dim[1];
  dim = INTEGER(GET_DIM(vmatrix)); nevent = dim[1];
  ntimes = LENGTH(times);

  ndeps = LENGTH(deps);
  if (ndeps > 0) didx = INTEGER(deps);

  PROTECT(Snames = GET_ROWNAMES(GET_DIMNAMES(xstart))); nprotect++;
  PROTECT(Pnames = GET_ROWNAMES(GET_DIMNAMES(params))); nprotect++;
  PROTECT(Cnames = GET_COLNAMES(GET_DIMNAMES(covar))); nprotect++;

  PROTECT(statenames = GET_SLOT(func,install("statenames"))); nprotect++;
  PROTECT(paramnames = GET_SLOT(func,install("paramnames"))); nprotect++;
  PROTECT(covarnames = GET_SLOT(func,install("covarnames"))); nprotect++;

  nstates = LENGTH(statenames);
  nparams = LENGTH(paramnames);
  ncovars = LENGTH(covarnames);

  PROTECT(fn = pomp_fun_handler(func,gnsi,&use_native)); nprotect++;

  if (use_native) {
    *((void **) (&(RXR))) = R_ExternalPtrAddr(fn);
  } else {
    RXR = (pomp_ssa_rate_fn *) default_ssa_rate_fn;
    PROTECT(RHO = (CLOENV(fn))); nprotect++;
    NVAR = nvar;
    NPAR = npar;
    PROTECT(JRATE = NEW_INTEGER(1)); nprotect++; // for internal use
    PROTECT(TIME = NEW_NUMERIC(1)); nprotect++;	// for internal use
    PROTECT(XVEC = NEW_NUMERIC(nvar)); nprotect++; // for internal use
    PROTECT(PVEC = NEW_NUMERIC(npar)); nprotect++; // for internal use
    PROTECT(CVEC = NEW_NUMERIC(covdim)); nprotect++; // for internal use
    SET_NAMES(XVEC,Snames); // make sure the names attribute is copied
    SET_NAMES(PVEC,Pnames); // make sure the names attribute is copied
    SET_NAMES(CVEC,Cnames); // make sure the names attribute is copied
    // set up the function call
    PROTECT(FCALL = LCONS(CVEC,args)); nprotect++;
    SET_TAG(FCALL,install("covars"));
    PROTECT(FCALL = LCONS(PVEC,FCALL)); nprotect++;
    SET_TAG(FCALL,install("params"));
    PROTECT(FCALL = LCONS(TIME,FCALL)); nprotect++;
    SET_TAG(FCALL,install("t"));
    PROTECT(FCALL = LCONS(XVEC,FCALL)); nprotect++;
    SET_TAG(FCALL,install("x"));
    PROTECT(FCALL = LCONS(JRATE,FCALL)); nprotect++;
    SET_TAG(FCALL,install("j"));
    PROTECT(FCALL = LCONS(fn,FCALL)); nprotect++;
    FIRST = 1;
  }

  xdim[0] = nvar; xdim[1] = nrep; xdim[2] = ntimes;
  PROTECT(X = makearray(3,xdim)); nprotect++;
  setrownames(X,Snames,3);

  if (nstates>0) {
    PROTECT(sindex = MATCHROWNAMES(xstart,statenames,"state variables")); nprotect++;
    sidx = INTEGER(sindex);
  } else {
    sidx = 0;
  }
  if (nparams>0) {
    PROTECT(pindex = MATCHROWNAMES(params,paramnames,"parameters")); nprotect++;
    pidx = INTEGER(pindex);
  } else {
    pidx = 0;
  }
  if (ncovars>0) {
    PROTECT(cindex = MATCHCOLNAMES(covar,covarnames,"covariates")); nprotect++;
    cidx = INTEGER(cindex);
  } else {
    cidx = 0;
  }
  if (nzeros>0) {
    PROTECT(zindex = MATCHROWNAMES(xstart,zeronames,"state variables")); nprotect++;
    zidx = INTEGER(zindex);
  } else {
    zidx = 0;
  }

  if (use_native) {
    set_pomp_userdata(args);
  }

  GetRNGstate();
  {
    int i;
    for (i = 0; i < nrep; i++) {
      SSA(RXR,i,nvar,nevent,npar,nrep,ntimes,
	  method,REAL(xstart),REAL(times),REAL(params),
	  REAL(X),REAL(e),REAL(vmatrix),REAL(dmatrix),ndeps,didx,
	  nzeros,zidx,sidx,pidx,ncovars,cidx,covlen,covdim,
	  REAL(tcovar),REAL(covar));
    }
  }
  PutRNGstate();

  if (use_native) {
    unset_pomp_userdata();
  }

  UNPROTECT(nprotect);
  return X;
}
