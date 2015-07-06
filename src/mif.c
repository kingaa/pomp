// -*- C++ -*-

#include "pomp_internal.h"
#include <Rdefines.h>

// implements the Ionides et al. (2006) MIF update rule

SEXP mif_update (SEXP pfp, SEXP theta, SEXP gamma, SEXP varfactor, 
		 SEXP sigma, SEXP pars)
{
  int nprotect = 0;
  double *v, *m1, *m2;
  double scal, sig, grad;
  int npar, ntimes, nfm, npv;
  SEXP FM, PV, newtheta;
  int *sidx, *thidx, *midx, *vidx, *dim;
  int i, j;

  sig = *(REAL(varfactor));
  scal = *(REAL(gamma))*(1+sig*sig);

  PROTECT(FM = GET_SLOT(pfp,install("filter.mean"))); nprotect++;
  PROTECT(PV = GET_SLOT(pfp,install("pred.var"))); nprotect++;

  npar = LENGTH(pars);
  dim = INTEGER(GET_DIM(FM)); nfm = dim[0]; ntimes = dim[1];
  dim = INTEGER(GET_DIM(PV)); npv = dim[0];

  sidx = INTEGER(PROTECT(MATCHNAMES(sigma,pars,"random-walk SDs"))); nprotect++;
  thidx = INTEGER(PROTECT(MATCHNAMES(theta,pars,"parameters"))); nprotect++;
  midx = INTEGER(PROTECT(MATCHROWNAMES(FM,pars,"filter-mean variables"))); nprotect++;
  vidx = INTEGER(PROTECT(MATCHROWNAMES(PV,pars,"prediction-variance variables"))); nprotect++;

  PROTECT(newtheta = duplicate(theta)); nprotect++;

  for (i = 0; i < npar; i++) {
    sig = REAL(sigma)[sidx[i]];
    m1 = REAL(theta)+thidx[i];
    m2 = REAL(FM)+midx[i];
    v = REAL(PV)+vidx[i];
    grad = (*m2-*m1)/(*v);
    for (j = 1, m1 = m2, m2 += nfm, v += npv; j < ntimes; j++, m1 = m2, m2 += nfm, v += npv) 
      grad += (*m2-*m1)/(*v);
    REAL(newtheta)[thidx[i]] += scal*sig*sig*grad;
  }

  UNPROTECT(nprotect);
  return newtheta;
}
