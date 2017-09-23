#include "pomp_internal.h"
#include <R_ext/Constants.h>

static void gillespie (double *t, double tmax, double *f, double *y, 
		       const double *v, int nvar, int nevent, Rboolean hasvname,
		       const int *ivmat) {
  double tstep, p;
  double vv;
  int i, j;

  // Determine time interval and update time
  double fsum = 0;
  for (j = 0; j < nevent; j++) fsum += f[j];

  if (fsum > 0.0) {
    tstep = exp_rand()/fsum;
    *t = *t+tstep;
  } else {
    *t = tmax;
    return;
  }

  if (*t >= tmax) {

    *t = tmax;

  } else {

    // Determine event, update pops & events
    p = fsum*unif_rand();
    int jevent = nevent-1;
    for (j = 0; j < jevent; j++) {
      if (p > f[j]) {
	p -= f[j];
      } else {
	jevent = j;
	break;
      }
    }

    for (i = 0; i < nvar; i++) {
      vv = v[i+nvar*jevent];
      if (vv != 0) {
	if (hasvname) {
	  y[ivmat[i]] += vv;
	} else {
	  y[i] += vv;
	}
      }
    }

  }
}

static void SSA (pomp_ssa_rate_fn *ratefun, int irep,
		 int nvar, int nevent, int npar, int nrep, int ntimes,
		 double *xstart, const double *times, const double *params, double *xout,
		 const double *v, int nzero, const int *izero,
		 const int *istate, const int *ipar, int ncovar, const int *icovar,
		 Rboolean hasvname, const int *ivmat,
		 int lcov, int mcov, double *tcov, double *cov, const double *hmax) {
  double t = times[0];
  double tmax;
  double *f = NULL;
  double *covars = NULL;
  const double *par;
  double y[nvar];
  struct lookup_table tab = {lcov, mcov, 0, tcov, cov};
  int i, j;

  if (mcov > 0) covars = (double *) Calloc(mcov,double);
  if (nevent > 0) f = (double *) Calloc(nevent,double);

  par = params+npar*irep;
  // Copy state variables
  memcpy(y,xstart+nvar*irep,nvar*sizeof(double));
  memcpy(xout+nvar*irep,y,nvar*sizeof(double));
  // Set appropriate states to zero
  for (i = 0; i < nzero; i++) y[izero[i]] = 0.0;
  // Initialize the covariate vector
  if (mcov > 0) table_lookup(&tab,t,covars);
  // Initialise propensity functions & tree
  for (j = 0; j < nevent; j++) {
    f[j] = ratefun(j+1,t,y,par,istate,ipar,icovar,mcov,covars);
    if (f[j] < 0.0)
      errorcall(R_NilValue,"'rate.fun' returns a negative rate");
  }

  int icount = 1;
  while (icount < ntimes) {

    R_CheckUserInterrupt();
    tmax = t + *hmax;
    tmax = (tmax > times[icount]) ? times[icount] : tmax;
    gillespie(&t,tmax,f,y,v,nvar,nevent,hasvname,ivmat);

    if (mcov > 0) table_lookup(&tab,t,covars);

    for (j = 0; j < nevent; j++) {
      f[j] = ratefun(j+1,t,y,par,istate,ipar,icovar,mcov,covars);
      if (f[j] < 0.0)
	errorcall(R_NilValue,"'rate.fun' returns a negative rate");
    }

    // Record output at required time points
    if (t >= times[icount]) {
      memcpy(xout+nvar*(irep+nrep*icount),y,nvar*sizeof(double));
      // Set appropriate states to zero at time of last observation
      for (i = 0; i < nzero; i++) y[izero[i]] = 0;
      icount++;
    }

  }

  if (mcov > 0) Free(covars);
  if (nevent > 0) Free(f);

}

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
      errorcall(R_NilValue,"user 'rates' must return a single scalar rate.");
    }
    FIRST = 0;
  }

  rate = *(REAL(AS_NUMERIC(ans)));
  UNPROTECT(nprotect);
  return rate;
}

SEXP SSA_simulator (SEXP func, SEXP xstart, SEXP times, SEXP params,
		    SEXP vmatrix, SEXP tcovar, SEXP covar,
		    SEXP zeronames, SEXP hmax, SEXP args, SEXP gnsi)
{
  int nprotect = 0;
  int *dim, xdim[3];
  int nvar, nvarv, nevent, npar, nrep, ntimes;
  int covlen, covdim;
  SEXP statenames, paramnames, covarnames;
  int nstates, nparams, ncovars;
  int nzeros = LENGTH(zeronames);
  pompfunmode use_native = undef;
  SEXP X, pindex, sindex, cindex, zindex, vindex;
  int *sidx, *pidx, *cidx, *zidx, *vidx;
  SEXP fn, Snames, Pnames, Cnames, Vnames;
  Rboolean hasvnames;

  dim = INTEGER(GET_DIM(xstart)); nvar = dim[0]; nrep = dim[1];
  dim = INTEGER(GET_DIM(params)); npar = dim[0];
  dim = INTEGER(GET_DIM(covar)); covlen = dim[0]; covdim = dim[1];
  dim = INTEGER(GET_DIM(vmatrix)); nvarv = dim[0]; nevent = dim[1];
  if (nvarv != nvar) {
    errorcall(R_NilValue,"number of state variables must equal the number of rows in v.");
  }
  ntimes = LENGTH(times);

  PROTECT(Snames = GET_ROWNAMES(GET_DIMNAMES(xstart))); nprotect++;
  PROTECT(Pnames = GET_ROWNAMES(GET_DIMNAMES(params))); nprotect++;
  PROTECT(Cnames = GET_COLNAMES(GET_DIMNAMES(covar))); nprotect++;
  PROTECT(Vnames = GET_ROWNAMES(GET_DIMNAMES(vmatrix))); nprotect++;

  PROTECT(statenames = GET_SLOT(func,install("statenames"))); nprotect++;
  PROTECT(paramnames = GET_SLOT(func,install("paramnames"))); nprotect++;
  PROTECT(covarnames = GET_SLOT(func,install("covarnames"))); nprotect++;

  nstates = LENGTH(statenames);
  nparams = LENGTH(paramnames);
  ncovars = LENGTH(covarnames);
  hasvnames = !isNull(Vnames);

  PROTECT(hmax = AS_NUMERIC(hmax)); nprotect++;

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
  if (hasvnames) {
    PROTECT(vindex = MATCHROWNAMES(xstart,Vnames,"state variables")); nprotect++;
    vidx = INTEGER(vindex);
  } else {
    vidx = 0;
  }

  if (use_native) {
    set_pomp_userdata(args);
  }

  GetRNGstate();
  {
    int i;
    for (i = 0; i < nrep; i++) {
      SSA(RXR,i,nvar,nevent,npar,nrep,ntimes,
	  REAL(xstart),REAL(times),REAL(params),
	  REAL(X),REAL(vmatrix),
	  nzeros,zidx,sidx,pidx,ncovars,cidx,hasvnames,vidx,covlen,covdim,
	  REAL(tcovar),REAL(covar),REAL(hmax));
    }
  }
  PutRNGstate();

  if (use_native) {
    unset_pomp_userdata();
  }

  UNPROTECT(nprotect);
  return X;
}
