#include "pomp_internal.h"

int gillespie (pomp_ssa_rate_fn *ratefun, double *t, double *f,
               double *y, const double *v, const double *d, const double *par,
               int nvar, int nevent, int npar,
               const int *istate, const int *ipar, int ncovar, const int *icovar,
               int mcov, const double *cov) {
  double tstep, p;
  int change[nvar];
  int i, j;
  // Determine time interval and update time
  double fsum = 0;
  for (j = 0; j < nevent; j++) fsum += f[j];
  p = unif_rand();
  if (fsum > 0.0) {
    tstep = -log(p)/fsum;
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
    if (v[i+nvar*jevent] != 0) {
      y[i] += v[i+nvar*jevent];
      change[i] = 1;
    }
  }
  // only updating events & tree entries that have changed
  for (j = 0; j < nevent; j++) {
    for (i = 0; i < nvar; i++) {
      if ((change[i] != 0) && (d[i*nvar*j] != 0)) {
        f[j] = (*ratefun)(j+1,*t,y,par,istate,ipar,icovar,mcov,cov);
	if (f[j] < 0.0) 
	  errorcall(R_NilValue,"'rate.fun' returns a negative rate");
        break;
      }
    }
  }
  return 0;
}

int kleap (pomp_ssa_rate_fn *ratefun, double kappa, double *t, double *f,
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
      if ((change[i] != 0) && (d[i*nvar*j] != 0)) {
        f[j] = (*ratefun)(j+1,*t,y,par,istate,ipar,icovar,mcov,cov);
	if (f[j] < 0.0) 
	  errorcall(R_NilValue,"'rate.fun' returns a negative rate");
        break;
      }
    }
  }
  return 0;
}

void SSA (pomp_ssa_rate_fn *ratefun, int irep,
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
  double covars[mcov], par[npar], y[nvar], f[nevent];
  struct lookup_table tab = {lcov, mcov, 0, tcov, cov};
  int i, j;

  // Copy parameters and states
  for (i = 0; i < npar; i++) par[i] = params[i+npar*irep];
  for (i = 0; i < nvar; i++)
    xout[i+irep*nvar] = y[i] = xstart[i+nvar*irep];
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
    if (method == 0) {	// Gillespie's next reaction method
      flag = gillespie(ratefun,&t,f,y,v,d,par,nvar,nevent,npar,istate,ipar,ncovar,icovar,mcov,covars);
    } else {	 // Cai's K-leap method
      // Determine kappa (most accurate but slowest method)
      double kappa, tmp;
      int k;
      for (i = 0, kappa = 1e9; i < ndeps; i++) {
        k = ideps[i];
        tmp = e[k]*y[k];
        kappa = (tmp < kappa) ? tmp : kappa;
        if (kappa < 2.0) break;
      }
      if (kappa < 2.0) {
        flag = gillespie(ratefun,&t,f,y,v,d,par,nvar,nevent,npar,istate,ipar,ncovar,icovar,mcov,covars);
      } else {
        kappa = floor(kappa);
        flag = kleap(ratefun,kappa,&t,f,y,v,d,par,nvar,nevent,npar,istate,ipar,ncovar,icovar,mcov,covars);
      }
    }
    // Record output at required time points
    while ((icount < ntimes) && (t >= times[icount])) {
      for (i = 0; i < nvar; i++)
        xout[i+nvar*(irep+nrep*icount)] = y[i];
      // Set appropriate states to zero
      for (i = 0; i < nzero; i++) y[izero[i]] =0.0;
      // Recompute if zero event-rate encountered
      if (flag) t = times[icount];
      icount++;
    }

    if ((mcov > 0) && (t <= tmax)) table_lookup(&tab,t,covars);
  }
}
