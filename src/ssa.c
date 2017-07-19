#include "pomp_internal.h"
#include <R_ext/Constants.h>

static int gillespie (pomp_ssa_rate_fn *ratefun, double *t, double *f,
		      double *y, const double *v, const double *d, const double *par,
		      int nvar, int nevent, int npar,
		      const int *istate, const int *ipar, int ncovar, const int *icovar,
		      int mcov, const double *cov) {
  double tstep, p;
  int change[nvar];
  double vv;
  int i, j;
  // Determine time interval and update time
  double fsum = 0;
  for (j = 0; j < nevent; j++) fsum += f[j];
  if (fsum > 0.0) {
    tstep = exp_rand()/fsum;
    *t = *t+tstep;
  } else {
    *t = R_PosInf;
    return 1;
  }
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
    change[i] = 0;
    vv = v[i+nvar*jevent];
    if (vv != 0) {
      y[i] += vv;
      change[i] = 1;
    }
  }
  // only updating events & tree entries that have changed
  for (j = 0; j < nevent; j++) {
    for (i = 0; i < nvar; i++) {
      if ((change[i] != 0) && (d[i+nvar*j] != 0)) {
        f[j] = (*ratefun)(j+1,*t,y,par,istate,ipar,icovar,mcov,cov);
        if (f[j] < 0.0)
          errorcall(R_NilValue,"'rate.fun' returns a negative rate");
        break;
      }
    }
  }
  return 0;
}

static int kleap (pomp_ssa_rate_fn *ratefun, double kappa, double *t, double *f,
		  double *y, const double *v, const double *d, const double *par,
		  int nvar, int nevent, int npar,
		  const int *istate, const int *ipar, int ncovar, const int *icovar,
		  int mcov, const double *cov) {
  double prob[nevent];
  int k[nevent];
  double kk, tstep;
  int change[nvar];
  int i, j;
  // Determine time interval and update time
  double fsum = 0;
  for (j = 0; j < nevent; j++) fsum += f[j];
  if (fsum > 0.0) {
    tstep = rgamma(kappa,1.0/fsum);
    *t = *t+tstep;
  } else {
    *t = R_PosInf;
    return 1;
  }
  // Determine frequency of events, update pops & events
  for (j = 0; j < nevent; j++) prob[j] = f[j]/fsum;
  rmultinom((int)kappa,prob,nevent,k);
  // some matrix-vector multiplication but only where necessary
  for (i = 0; i < nvar; i++) change[i] = 0;
  for (j = 0; j < nevent; j++) {
    if (k[j] != 0) {
      kk = (double) k[j];
      for (i = 0; i < nvar; i++) {
        if (v[i+nvar*j] != 0) {
          y[i] += kk*v[i+nvar*j];
          change[i] = 1;
        }
      }
    }
  }
  // only updating events & tree entries that have changed
  for (j = 0; j < nevent; j++) {
    for (i = 0; i < nvar; i++) {
      if ((change[i] != 0) && (d[i+nvar*j] != 0)) {
        f[j] = (*ratefun)(j+1,*t,y,par,istate,ipar,icovar,mcov,cov);
        if (f[j] < 0.0)
          errorcall(R_NilValue,"'rate.fun' returns a negative rate");
        break;
      }
    }
  }
  return 0;
}

static void SSA (pomp_ssa_rate_fn *ratefun, int irep,
		 int nvar, int nevent, int npar, int nrep, int ntimes,
		 int method,
		 double *xstart, const double *times, const double *params, double *xout,
		 const double *e, const double *v, const double *d,
		 int ndeps, const int *ideps, int nzero, const int *izero,
		 const int *istate, const int *ipar, int ncovar, const int *icovar,
		 int lcov, int mcov, double *tcov, double *cov) {
  int flag = 0;
  double t = times[0];
  double tmax = times[ntimes-1];
  double *covars = NULL;
  double *f = NULL;
  double par[npar], y[nvar], ynext[nvar];
  struct lookup_table tab = {lcov, mcov, 0, tcov, cov};
  int i, j;

  if (mcov > 0) covars = (double *) Calloc(mcov,double);
  if (nevent > 0) f = (double *) Calloc(nevent,double);

  // Copy parameters and states
  memcpy(par,params+npar*irep,npar*sizeof(double));
  memcpy(y,xstart+nvar*irep,nvar*sizeof(double));
  memcpy(xout+nvar*irep,xstart+nvar*irep,nvar*sizeof(double));
  // Set appropriate states to zero
  for (i = 0; i < nzero; i++) y[izero[i]] = 0.0;
  memcpy(ynext,y,nvar*sizeof(double));
  // Initialize the covariate vector
  if (mcov > 0) table_lookup(&tab,t,covars);
  // Initialise propensity functions & tree
  for (j = 0; j < nevent; j++) {
    f[j] = ratefun(j+1,t,ynext,par,istate,ipar,icovar,mcov,covars);
    if (f[j] < 0.0)
      errorcall(R_NilValue,"'rate.fun' returns a negative rate");
  }
  int icount = 1;
  while (icount < ntimes) {
    R_CheckUserInterrupt();
    if (method == 0) {	// Gillespie's next reaction method
      flag = gillespie(ratefun,&t,f,ynext,v,d,par,nvar,nevent,npar,istate,ipar,
        ncovar,icovar,mcov,covars);
    } else {	 // Cai's K-leap method
      // Determine kappa (most accurate but slowest method)
      double kappa, tmp;
      int k;
      for (i = 0, kappa = 1e9; i < ndeps; i++) {
        k = ideps[i];
        tmp = e[k]*ynext[k];
        kappa = (tmp < kappa) ? tmp : kappa;
        if (kappa < 2.0) break;
      }
      if (kappa < 2.0) {
        flag = gillespie(ratefun,&t,f,ynext,v,d,par,nvar,nevent,npar,istate,
          ipar,ncovar,icovar,mcov,covars);
      } else {
        kappa = floor(kappa);
        flag = kleap(ratefun,kappa,&t,f,ynext,v,d,par,nvar,nevent,npar,istate,
          ipar,ncovar,icovar,mcov,covars);
      }
    }

    // Record output at required time points
    if (t >= times[icount]) {
      while ((icount < ntimes) && (t >= times[icount])) {
        memcpy(xout+nvar*(irep+nrep*icount),y,nvar*sizeof(double));
        // Set appropriate states to zero at time of last observation
        for (i = 0; i < nzero; i++) y[izero[i]] = 0;
        // Recompute if zero event-rate encountered
        if (flag) t = times[icount];
        icount++;
      }
      memcpy(y,ynext,nvar*sizeof(double));
      for (i = 0; i < nzero; i++) ynext[izero[i]] = 0;
    }

    if ((mcov > 0) && (t <= tmax)) table_lookup(&tab,t,covars);

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
      UNPROTECT(nprotect);
      errorcall(R_NilValue,"user 'rates' must return a single scalar rate.");
    }
    FIRST = 0;
  }

  rate = *(REAL(AS_NUMERIC(ans)));
  UNPROTECT(nprotect);
  return rate;
}

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
  pompfunmode use_native = undef;
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
