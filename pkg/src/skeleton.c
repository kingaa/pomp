// dear emacs, please treat this as -*- C++ -*-

#include "pomp_internal.h"
#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

static void eval_skel (pomp_skeleton *vf,
		       double *f, 
		       double *x, double *times, double *params, 
		       int *ndim,
		       int *stateindex, int *parindex, int *covindex,
		       double *time_table, double *covar_table)
{
  double t, *xp, *pp, *fp;
  int nvar = ndim[0];
  int npar = ndim[1];
  int nrepp = ndim[2];
  int nrepx = ndim[3];
  int ntimes = ndim[4];
  int covlen = ndim[5];
  int covdim = ndim[6];
  double covar_fn[covdim];
  int nrep;
  int k, p;
  
  // set up the covariate table
  struct lookup_table covariate_table = {covlen, covdim, 0, time_table, covar_table};

  nrep = (nrepp > nrepx) ? nrepp : nrepx;
  
  for (k = 0; k < ntimes; k++) { // loop over times

    R_CheckUserInterrupt();	// check for user interrupt

    t = times[k];

    // interpolate the covar functions for the covariates
    if (covdim > 0) 
      table_lookup(&covariate_table,t,covar_fn,0);
    
    for (p = 0; p < nrep; p++) { // loop over replicates
      
      fp = &f[nvar*(p+nrep*k)];
      xp = &x[nvar*((p%nrepx)+nrepx*k)];
      pp = &params[npar*(p%nrepp)];

      (*vf)(fp,xp,pp,stateindex,parindex,covindex,covdim,covar_fn,t);
      
    }
  }
}

// these global objects will pass the needed information to the user-defined function (see 'default_skel_fn')
// each of these is allocated once, globally, and refilled many times
static SEXP _pomp_skel_Xvec;	// state variable vector
static SEXP _pomp_skel_Pvec;	// parameter vector
static SEXP _pomp_skel_Cvec;	// covariate vector
static SEXP _pomp_skel_time;	// time
static int  _pomp_skel_nvar;	// number of state variables
static int  _pomp_skel_npar;	// number of parameters
static SEXP _pomp_skel_envir;	// function's environment
static SEXP _pomp_skel_fcall;	// function call
static SEXP _pomp_skel_vnames;	// names of state variables
static int *_pomp_skel_vindex;	// indices of state variables
static int  _pomp_skel_first;	// first evaluation?

#define XVEC    (_pomp_skel_Xvec)
#define PVEC    (_pomp_skel_Pvec)
#define CVEC    (_pomp_skel_Cvec)
#define TIME    (_pomp_skel_time)
#define NVAR    (_pomp_skel_nvar)
#define NPAR    (_pomp_skel_npar)
#define RHO     (_pomp_skel_envir)
#define FCALL   (_pomp_skel_fcall)
#define VNAMES  (_pomp_skel_vnames)
#define VIDX    (_pomp_skel_vindex)
#define FIRST   (_pomp_skel_first)

// this is the vectorfield that is evaluated when the user supplies an R function
// (and not a native routine)
// Note that stateindex, parindex, covindex are ignored.
static void default_skel_fn (double *f, double *x, double *p, 
			     int *stateindex, int *parindex, int *covindex, 
			     int covdim, double *covar, double t)
{
  int nprotect = 0;
  int k;
  double *xp;
  SEXP ans, nm, idx;
  int use_names = 0, *op;

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
    if (LENGTH(ans)!=NVAR)
      error("user 'skeleton' returns a vector of %d state variables but %d are expected",LENGTH(ans),NVAR);
    PROTECT(nm = GET_NAMES(ans)); nprotect++;
    use_names = !isNull(nm);
    if (use_names) {	  // match names against known names of states
      PROTECT(idx = matchnames(VNAMES,nm)); nprotect++;
      op = INTEGER(idx);
      for (k = 0; k < NVAR; k++) VIDX[k] = op[k];
    } else {
      VIDX = 0;
    }
    FIRST = 0;
  }

  xp = REAL(AS_NUMERIC(ans));
  if (VIDX == 0) {
    for (k = 0; k < NVAR; k++) f[k] = xp[k];
  } else {
    for (k = 0; k < NVAR; k++) f[VIDX[k]] = xp[k];
  }

  UNPROTECT(nprotect);
}

SEXP do_skeleton (SEXP object, SEXP x, SEXP t, SEXP params, SEXP fun)
{
  int nprotect = 0;
  int *dim, nvars, npars, nrepsp, nrepsx, nreps, ntimes, covlen, covdim;
  int use_native;
  int nstates, nparams, ncovars;
  double *xp;
  SEXP dimP, dimX, fn, F;
  SEXP tcovar, covar;
  SEXP statenames, paramnames, covarnames;
  SEXP sindex, pindex, cindex;
  SEXP Pnames, Cnames;
  pomp_skeleton *ff = NULL;
  int k, len;

  PROTECT(t = AS_NUMERIC(t)); nprotect++;
  ntimes = LENGTH(t);

  PROTECT(x = as_state_array(x)); nprotect++;
  dim = INTEGER(GET_DIM(x));
  nvars = dim[0]; nrepsx = dim[1];
  if (ntimes != dim[2])
    error("skeleton error: length of 't' and 3rd dimension of 'x' do not agree");

  PROTECT(params = as_matrix(params)); nprotect++;
  dim = INTEGER(GET_DIM(params));
  npars = dim[0]; nrepsp = dim[1];

  nreps = (nrepsp > nrepsx) ? nrepsp : nrepsx;

  if ((nreps % nrepsp != 0) || (nreps % nrepsx != 0))
    error("skeleton error: larger number of replicates is not a multiple of smaller");

  PROTECT(tcovar =  GET_SLOT(object,install("tcovar"))); nprotect++;
  PROTECT(covar =  GET_SLOT(object,install("covar"))); nprotect++;

  PROTECT(statenames =  GET_SLOT(object,install("statenames"))); nprotect++;
  PROTECT(paramnames =  GET_SLOT(object,install("paramnames"))); nprotect++;
  PROTECT(covarnames =  GET_SLOT(object,install("covarnames"))); nprotect++;
  nstates = LENGTH(statenames);
  nparams = LENGTH(paramnames);
  ncovars = LENGTH(covarnames);

  dim = INTEGER(GET_DIM(covar)); 
  covlen = dim[0]; covdim = dim[1];
  PROTECT(VNAMES = GET_ROWNAMES(GET_DIMNAMES(x))); nprotect++;
  PROTECT(Pnames = GET_ROWNAMES(GET_DIMNAMES(params))); nprotect++;
  PROTECT(Cnames = GET_COLNAMES(GET_DIMNAMES(covar))); nprotect++;

  PROTECT(fn = VECTOR_ELT(fun,0)); nprotect++;
  use_native = INTEGER(VECTOR_ELT(fun,1))[0];
  
  switch (use_native) {
  case 0:			// R skeleton
    ff = (pomp_skeleton *) default_skel_fn;
    PROTECT(RHO = (CLOENV(fn))); nprotect++;
    NVAR = nvars;			// for internal use
    NPAR = npars;			// for internal use
    PROTECT(TIME = NEW_NUMERIC(1)); nprotect++;	// for internal use
    PROTECT(XVEC = NEW_NUMERIC(nvars)); nprotect++; // for internal use
    xp = REAL(XVEC);
    for (k = 0; k < nvars; k++) xp[k] = 0.0;
    PROTECT(PVEC = NEW_NUMERIC(npars)); nprotect++; // for internal use
    PROTECT(CVEC = NEW_NUMERIC(covdim)); nprotect++; // for internal use
    SET_NAMES(XVEC,VNAMES); // make sure the names attribute is copied
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
    VIDX = (int *) R_alloc(nvars,sizeof(int));
    FIRST = 1;
    break;
  case 1:			// native skeleton
    ff = (pomp_skeleton *) R_ExternalPtrAddr(fn);
    VIDX = 0;
    break;
  default:
    error("unrecognized 'use' slot in 'skeleton'");
    break;
  }

  {
    int fdim[3];
    fdim[0] = nvars; fdim[1] = nreps; fdim[2] = ntimes;
    PROTECT(F = makearray(3,fdim)); nprotect++; 
    setrownames(F,VNAMES,3);
    xp = REAL(F);
    //  for (k = 0, len = nvars*nreps*ntimes; k < len; k++) xp[k] = 0.0;
  }

  if (nstates > 0) {
    PROTECT(sindex = MATCHROWNAMES(x,statenames)); nprotect++;
  } else {
    PROTECT(sindex = NEW_INTEGER(0)); nprotect++;
  }
  if (nparams > 0) {
    PROTECT(pindex = MATCHROWNAMES(params,paramnames)); nprotect++;
  } else {
    PROTECT(pindex = NEW_INTEGER(0)); nprotect++;
  }
  if (ncovars > 0) {
    PROTECT(cindex = MATCHCOLNAMES(covar,covarnames)); nprotect++;
  } else {
    PROTECT(cindex = NEW_INTEGER(0)); nprotect++;
  }

  {
    int ndim[7];
    ndim[0] = nvars; ndim[1] = npars; ndim[2] = nrepsp; ndim[3] = nrepsx; 
    ndim[4] = ntimes; ndim[5] = covlen; ndim[6] = covdim;

    eval_skel(ff,REAL(F),REAL(x),REAL(t),REAL(params),
	      ndim,INTEGER(sindex),INTEGER(pindex),INTEGER(cindex),
	      REAL(tcovar),REAL(covar));
  }

  UNPROTECT(nprotect);
  return F;
}

#undef XVEC
#undef PVEC
#undef CVEC
#undef TIME
#undef NVAR
#undef NPAR  
#undef RHO   
#undef FCALL 
#undef VNAMES


static SEXP *_pomp_vf_eval_object;
static SEXP *_pomp_vf_eval_params;
static SEXP *_pomp_vf_eval_skelfun;
static SEXP *_pomp_vf_eval_xnames;
static int _pomp_vf_eval_xdim[3];
#define OBJECT        (_pomp_vf_eval_object)
#define PARAMS        (_pomp_vf_eval_params)
#define SKELFUN       (_pomp_vf_eval_skelfun)
#define XNAMES        (_pomp_vf_eval_xnames)
#define XDIM          (_pomp_vf_eval_xdim)

void pomp_desolve_init (SEXP object, SEXP params, SEXP fun, SEXP statenames, SEXP nvar, SEXP nrep) {
  OBJECT = &object;
  PARAMS = &params;
  SKELFUN = &fun;
  XNAMES = &statenames;
  XDIM[0] = INTEGER(AS_INTEGER(nvar))[0];
  XDIM[1] = INTEGER(AS_INTEGER(nrep))[0];
  XDIM[2] = 1;
}


void pomp_vf_eval (int *neq, double *t, double *y, double *ydot, double *yout, int *ip) 
{
  SEXP T, X, dXdt;
  int dim[3];
  
  PROTECT(T = NEW_NUMERIC(1));
  PROTECT(X = makearray(3,XDIM));
  setrownames(X,*(XNAMES),3);
  REAL(T)[0] = *t;
  memcpy(REAL(X),y,(*neq)*sizeof(double));
  PROTECT(dXdt = do_skeleton(*(OBJECT),X,T,*(PARAMS),*(SKELFUN)));
  memcpy(ydot,REAL(dXdt),(*neq)*sizeof(double));
  
  UNPROTECT(3);
}

#undef XDIM
#undef XNAMES
#undef SKELFUN
#undef PARAMS
#undef OBJECT
