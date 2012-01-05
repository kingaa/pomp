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

SEXP iterate_map (SEXP object, SEXP times, SEXP t0, SEXP x0, SEXP params,
		  SEXP dt, SEXP zeronames)
{
  int nprotect = 0;
  SEXP ans;
  SEXP x, f, skel, time, zindex;
  int nvars, nreps, ntimes;
  int nzeros, nsteps;
  int h, i, j, k;
  int ndim[3];
  double *tm, *tp, *xp, *fp, *ap;
  int *zidx;
  double deltat;

  PROTECT(t0 = AS_NUMERIC(t0)); nprotect++;
  PROTECT(x = duplicate(AS_NUMERIC(x0))); nprotect++;
  xp = REAL(x);

  PROTECT(times = AS_NUMERIC(times)); nprotect++;
  ntimes = LENGTH(times);
  tp = REAL(times);

  nvars = INTEGER(GET_DIM(x0))[0];
  nreps = INTEGER(GET_DIM(x0))[1];
  if (nreps != INTEGER(GET_DIM(params))[1])
    error("mismatch in dimensions of 'x0' and 'params'");

  PROTECT(time = NEW_NUMERIC(1)); nprotect++;
  tm = REAL(time);
  *tm = REAL(t0)[0];

  ndim[0] = nvars; ndim[1] = nreps; ndim[2] = ntimes;
  PROTECT(ans = makearray(3,ndim)); nprotect++;
  setrownames(ans,GET_ROWNAMES(GET_DIMNAMES(x0)),3);
  ap = REAL(ans);

  PROTECT(skel = get_pomp_fun(GET_SLOT(object,install("skeleton")))); nprotect++;

  PROTECT(dt = AS_NUMERIC(dt)); nprotect++;
  if (REAL(dt)[0]<=0) 
    error("timestep 'dt' must be positive!");

  PROTECT(zeronames = AS_CHARACTER(zeronames)); nprotect++;
  nzeros = LENGTH(zeronames);

  if (nzeros>0) {
    PROTECT(zindex = MATCHROWNAMES(x0,zeronames)); nprotect++;
    zidx = INTEGER(zindex);
  } else {
    zidx = 0;
  }

  if (nzeros>0) {
    for (j = 0; j < nreps; j++) {
      for (i = 0; i < nzeros; i++) {
	xp[zidx[i]+nvars*j] = 0.0;
      }
    }
  }

  for (k = 0; k < ntimes; k++, tp++) {

    if (k > 0) *tm = *(tp-1);
    // compute number of steps to take
    deltat = REAL(dt)[0];
    nsteps = num_euler_steps(*tm,*tp,&deltat); 

    for (h = 0; h < nsteps; h++) {
      PROTECT(f = do_skeleton(object,x,time,params,skel));
      fp = REAL(f);
      for (j = 0; j < nreps; j++) {
	for (i = 0; i < nvars; i++) {
	  xp[i+nvars*j] = fp[i+nvars*j];
	}
      }
      UNPROTECT(1);
      *tm += deltat;
      if (h == nsteps-2) {	// penultimate step
	deltat = *(tp+1)-*tm;
	*tm = *(tp+1)-deltat;
      }
    }
    for (j = 0; j < nreps; j++) {
      for (i = 0; i < nvars; i++) {
	ap[i+nvars*(j+nreps*k)] = xp[i+nvars*j];
      }
      for (i = 0; i < nzeros; i++) {
	xp[zidx[i]+nvars*j] = 0.0;
      }
    }
  }

  UNPROTECT(nprotect);
  return ans;
}
