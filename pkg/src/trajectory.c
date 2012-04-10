// -*- C++ -*-

#include <Rdefines.h>
#include <R_ext/Constants.h>

#include "pomp_internal.h"

// copy t(x[-1,-1]) -> y[,rep,]
SEXP traj_transp_and_copy (SEXP y, SEXP x, SEXP rep) {
  int nprotect = 0;
  SEXP ans = R_NilValue;
  int nvars, nreps, ntimes;
  int i, j, k, m, n;
  double *xp, *yp;

  j = INTEGER(rep)[0]-1;
  nvars = INTEGER(GET_DIM(y))[0];
  nreps = INTEGER(GET_DIM(y))[1];
  ntimes = INTEGER(GET_DIM(y))[2];
  n = INTEGER(GET_DIM(x))[0];
  m = nvars*nreps;

  for (i = 0, xp = REAL(x)+n+1; i < nvars; i++, xp += n) {
    for (k = 0, yp = REAL(y)+i+nvars*j; k < ntimes; k++, yp += m) {
      *yp = xp[k];
    }
  }

  UNPROTECT(nprotect);
  return ans;
}

SEXP iterate_map (SEXP object, SEXP times, SEXP t0, SEXP x0, SEXP params)
{
  int nprotect = 0;
  SEXP ans;
  SEXP x, skel, time, zeronames;
  int nvars, nreps, ntimes;
  int nzeros;
  int nsteps;
  int h, i, j, k;
  int ndim[3];
  double *tm, *tp, *xp, *fp, *ap;
  int *zidx;
  double deltat;

  PROTECT(x = as_state_array(duplicate(AS_NUMERIC(x0)))); nprotect++;
  xp = REAL(x);

  PROTECT(times = AS_NUMERIC(times)); nprotect++;
  ntimes = LENGTH(times);
  tp = REAL(times);

  nvars = INTEGER(GET_DIM(x0))[0];
  nreps = INTEGER(GET_DIM(x0))[1];
  if (nreps != INTEGER(GET_DIM(params))[1])
    error("mismatch in dimensions of 'x0' and 'params'");

  PROTECT(time = duplicate(AS_NUMERIC(t0))); nprotect++;
  tm = REAL(time);

  ndim[0] = nvars; ndim[1] = nreps; ndim[2] = ntimes;
  PROTECT(ans = makearray(3,ndim)); nprotect++;
  setrownames(ans,GET_ROWNAMES(GET_DIMNAMES(x0)),3);
  ap = REAL(ans);

  PROTECT(skel = get_pomp_fun(GET_SLOT(object,install("skeleton")))); nprotect++;
  deltat = *(REAL(GET_SLOT(object,install("skelmap.delta.t"))));

  PROTECT(zeronames = GET_SLOT(object,install("zeronames"))); nprotect++;
  nzeros = LENGTH(zeronames);
  if (nzeros>0) {
    zidx = INTEGER(PROTECT(MATCHROWNAMES(x0,zeronames))); nprotect++;
  } else {
    zidx = 0;
  }

  if (nzeros>0) {
    for (j = 0; j < nreps; j++)
      for (i = 0; i < nzeros; i++) xp[zidx[i]+nvars*j] = 0.0;
  }

  for (k = 0; k < ntimes; k++, tp++) {

    nsteps = num_map_steps(*tm,*tp,deltat); 

    for (h = 0; h < nsteps; h++) {
      fp = REAL(do_skeleton(object,x,time,params,skel));
      for (j = 0; j < nreps; j++)
	for (i = 0; i < nvars; i++) xp[i+nvars*j] = fp[i+nvars*j];
      *tm += deltat;
    }

    for (j = 0; j < nreps; j++) {
      for (i = 0; i < nvars; i++) ap[i+nvars*(j+nreps*k)] = xp[i+nvars*j];
      for (i = 0; i < nzeros; i++) xp[zidx[i]+nvars*j] = 0.0;
    }

  }

  UNPROTECT(nprotect);
  return ans;
}

static struct {
  SEXP object;
  SEXP params;
  SEXP skelfun;
  SEXP xnames;
  int xdim[3];
} _pomp_vf_eval_common;


#define COMMON(X)    (_pomp_vf_eval_common.X)

void pomp_desolve_init (SEXP object, SEXP params, SEXP fun, SEXP statenames, SEXP nvar, SEXP nrep) {
  COMMON(object) = object;
  COMMON(params) = params;
  COMMON(skelfun) = fun;
  COMMON(xnames) = statenames;
  COMMON(xdim)[0] = INTEGER(AS_INTEGER(nvar))[0];
  COMMON(xdim)[1] = INTEGER(AS_INTEGER(nrep))[0];
  COMMON(xdim)[2] = 1;
}


void pomp_vf_eval (int *neq, double *t, double *y, double *ydot, double *yout, int *ip) 
{
  SEXP T, X, dXdt;
  
  PROTECT(T = NEW_NUMERIC(1));
  PROTECT(X = makearray(3,COMMON(xdim)));
  setrownames(X,COMMON(xnames),3);
  REAL(T)[0] = *t;
  memcpy(REAL(X),y,(*neq)*sizeof(double));
  PROTECT(dXdt = do_skeleton(COMMON(object),X,T,COMMON(params),COMMON(skelfun)));
  memcpy(ydot,REAL(dXdt),(*neq)*sizeof(double));
  
  UNPROTECT(3);
}

#undef COMMON
