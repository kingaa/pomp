// dear emacs, please treat this as -*- C++ -*-

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Constants.h>
#include <string.h>

#include "pomp_internal.h"


static R_INLINE SEXP add_args (SEXP args, SEXP Snames, SEXP Pnames, SEXP Cnames)
{
  int nprotect = 0;
  SEXP var;
  int v;

  // we construct the call from end to beginning
  // covariates, parameter, states, then time and 'j'

  // Covariates
  for (v = LENGTH(Cnames)-1; v >= 0; v--) {
    PROTECT(var = NEW_NUMERIC(1)); nprotect++;
    PROTECT(args = LCONS(var,args)); nprotect++;
    SET_TAG(args,install(CHAR(STRING_ELT(Cnames,v))));
  }

  // Parameters
  for (v = LENGTH(Pnames)-1; v >= 0; v--) {
    PROTECT(var = NEW_NUMERIC(1)); nprotect++;
    PROTECT(args = LCONS(var,args)); nprotect++;
    SET_TAG(args,install(CHAR(STRING_ELT(Pnames,v))));
  }

  // Latent state variables
  for (v = LENGTH(Snames)-1; v >= 0; v--) {
    PROTECT(var = NEW_NUMERIC(1)); nprotect++;
    PROTECT(args = LCONS(var,args)); nprotect++;
    SET_TAG(args,install(CHAR(STRING_ELT(Snames,v))));
  }

  // Time
  PROTECT(var = NEW_NUMERIC(1)); nprotect++;
  PROTECT(args = LCONS(var,args)); nprotect++;
  SET_TAG(args,install("t"));

  // 'j'
  PROTECT(var = NEW_INTEGER(1)); nprotect++;
  PROTECT(args = LCONS(var,args)); nprotect++;
  SET_TAG(args,install("j"));

  UNPROTECT(nprotect);
  return args;

}

static void gillespie (double *t, double tmax, double *f, double *y,
  const double *v, int nvar, int nevent, Rboolean hasvname,
  const int *ivmat)
{
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
  int mcov, lookup_table_t *tab, const double *hmax)
{
  double t = times[0];
  double tmax;
  double *f = NULL;
  double *covars = NULL;
  const double *par;
  double y[nvar];
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
  if (mcov > 0) table_lookup(tab,t,covars);
  // Initialise propensity functions & tree
  for (j = 0; j < nevent; j++) {
    f[j] = ratefun(j+1,t,y,par,istate,ipar,icovar,covars);
    if (f[j] < 0.0)
      errorcall(R_NilValue,"'rate.fun' returns a negative rate");
  }

  int icount = 1;
  while (icount < ntimes) {

    R_CheckUserInterrupt();
    tmax = t + *hmax;
    tmax = (tmax > times[icount]) ? times[icount] : tmax;
    gillespie(&t,tmax,f,y,v,nvar,nevent,hasvname,ivmat);

    if (mcov > 0) table_lookup(tab,t,covars);

    for (j = 0; j < nevent; j++) {
      f[j] = ratefun(j+1,t,y,par,istate,ipar,icovar,covars);
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

// these global objects will pass the needed information to the user-defined function (see '__pomp_Rfun_ssa_ratefn')
// each of these is allocated once, globally, and refilled many times
static int  __ssa_nvar;	// number of state variables
static int  __ssa_npar;	// number of parameters
static int  __ssa_ncov;	// number of covariates
static SEXP __ssa_args;	// function arguments
static SEXP __ssa_ratefn;	// function itself
static int  __ssa_first;	// first evaluation?
static pomp_ssa_rate_fn *__ssa_rxrate; // function computing reaction rates

#define NVAR    (__ssa_nvar)
#define NPAR    (__ssa_npar)
#define NCOV    (__ssa_ncov)
#define ARGS    (__ssa_args)
#define RATEFN  (__ssa_ratefn)
#define FIRST   (__ssa_first)
#define RXR     (__ssa_rxrate)

static double __pomp_Rfun_ssa_ratefn (int j, double t, const double *x, const double *p,
  int *stateindex, int *parindex, int *covindex, double *c)
{
  SEXP var = ARGS, ans;
  int v;
  double rate;

  *(INTEGER(CAR(var))) = j; var = CDR(var);
  *(REAL(CAR(var))) = t; var = CDR(var);
  for (v = 0; v < NVAR; v++, x++, var=CDR(var)) *(REAL(CAR(var))) = *x;
  for (v = 0; v < NPAR; v++, p++, var=CDR(var)) *(REAL(CAR(var))) = *p;
  for (v = 0; v < NCOV; v++, c++, var=CDR(var)) *(REAL(CAR(var))) = *c;

  PROTECT(ans = eval(LCONS(RATEFN,ARGS),CLOENV(RATEFN)));

  if (FIRST) {
    if (LENGTH(ans) != 1)
      errorcall(R_NilValue,"'rate.fun' must return a single numeric rate.");
    FIRST = 0;
  }

  rate = *(REAL(AS_NUMERIC(ans)));
  UNPROTECT(1);
  return rate;
}

SEXP SSA_simulator (SEXP func, SEXP xstart, SEXP times, SEXP params,
  SEXP vmatrix, SEXP covar, SEXP zeronames, SEXP hmax, SEXP args,
  SEXP gnsi)
{
  int nprotect = 0;
  int *dim, xdim[3];
  int nvar, nvarv, nevent, npar, nrep, ntimes;
  SEXP statenames, paramnames, covarnames;
  int nstates, nparams, ncovars, covdim;
  int nzeros = LENGTH(zeronames);
  pompfunmode mode = undef;
  SEXP X, pindex, sindex, cindex, zindex, vindex;
  int *sidx, *pidx, *cidx, *zidx, *vidx;
  SEXP fn, Snames, Pnames, Cnames, Vnames;
  lookup_table_t covariate_table;
  Rboolean hasvnames;

  dim = INTEGER(GET_DIM(xstart)); nvar = dim[0]; nrep = dim[1];
  dim = INTEGER(GET_DIM(params)); npar = dim[0];
  dim = INTEGER(GET_DIM(vmatrix)); nvarv = dim[0]; nevent = dim[1];

  if (nvarv != nvar)
    errorcall(R_NilValue,"number of state variables must equal the number of rows in v.");

  ntimes = LENGTH(times);

  PROTECT(Snames = GET_ROWNAMES(GET_DIMNAMES(xstart))); nprotect++;
  PROTECT(Pnames = GET_ROWNAMES(GET_DIMNAMES(params))); nprotect++;
  PROTECT(Cnames = get_covariate_names(covar)); nprotect++;
  PROTECT(Vnames = GET_ROWNAMES(GET_DIMNAMES(vmatrix))); nprotect++;

  covariate_table = make_covariate_table(covar,&covdim);

  PROTECT(statenames = GET_SLOT(func,install("statenames"))); nprotect++;
  PROTECT(paramnames = GET_SLOT(func,install("paramnames"))); nprotect++;
  PROTECT(covarnames = GET_SLOT(func,install("covarnames"))); nprotect++;

  nstates = LENGTH(statenames);
  nparams = LENGTH(paramnames);
  ncovars = LENGTH(covarnames);
  hasvnames = !isNull(Vnames);

  PROTECT(hmax = AS_NUMERIC(hmax)); nprotect++;

  PROTECT(fn = pomp_fun_handler(func,gnsi,&mode)); nprotect++;

  switch (mode) {

  case undef: default:

    errorcall(R_NilValue,"unrecognized 'mode' %d",mode); // # nocov

  case native: case regNative:

    *((void **) (&(RXR))) = R_ExternalPtrAddr(fn);

    set_pomp_userdata(args);

    break;

  case Rfun:

    RXR = (pomp_ssa_rate_fn *) __pomp_Rfun_ssa_ratefn;
    NVAR = nvar;
    NPAR = npar;
    NCOV = covdim;
    PROTECT(ARGS = add_args(args,Snames,Pnames,Cnames)); nprotect++;
    PROTECT(RATEFN = fn); nprotect++;
    FIRST = 1;

    break;

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
    PROTECT(cindex = matchnames(Cnames,covarnames,"covariates")); nprotect++;
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

  GetRNGstate();
  {
    int i;
    for (i = 0; i < nrep; i++) {
      SSA(RXR,i,nvar,nevent,npar,nrep,ntimes,
        REAL(xstart),REAL(times),REAL(params),
        REAL(X),REAL(vmatrix),
        nzeros,zidx,sidx,pidx,ncovars,cidx,hasvnames,vidx,
        covdim,&covariate_table,REAL(hmax));
    }
  }
  PutRNGstate();

  if (mode == native || mode == regNative) {
    unset_pomp_userdata();
  }

  UNPROTECT(nprotect);
  return X;
}
