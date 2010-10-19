// -*- mode: C++ -*-

#include "pomp_internal.h"
#include "pomp_mat.h"
#include <string.h>

// Campbell's robust variance-covariance estimation approach
// described on p. 231 of Krzanowski 1988
// with additional pre-conditioning for numerical stability
// translated into C from R code by Simon N. Wood

void robust_synth_loglik (double *y, int *dim, double *ydat, double *loglik) {
  // 'ydat' is destroyed
  // 'y' is preserved
  int nrow = dim[0];
  int ncol = dim[1];
  double alpha = 2.0, beta = 1.25;
  double x, xx, *yp, wbar;
  double w[nrow], tau[ncol];
  int pivot[nrow];
  double *y1, *y2;
  double d, d0, rss, half_log_det;
  int i, j, k;

  half_log_det = ncol*M_LN_SQRT_2PI;

  y1 = (double *) R_alloc(nrow*ncol,sizeof(double));
  y2 = (double *) R_alloc(nrow*ncol,sizeof(double));

  // compute column means, center each column, compute preconditioner
  memcpy(y1,y,nrow*ncol*sizeof(double));
  for (yp = y1, j = 0; j < ncol; j++, yp += nrow) {
    for (xx = 0, i = 0; i < nrow; i++) xx += yp[i];
    xx /= nrow;
    for (i = 0; i < nrow; i++) yp[i] -= xx; // center the column
    for (xx = 0, i = 0; i < nrow; i++) xx += yp[i]*yp[i];
    x = sqrt(xx/nrow);		   // column SD
    for (i = 0; i < nrow; i++) yp[i] /= x; // precondition
  }

  // do first QR decomposition & backsolve

  memcpy(y2,y1,nrow*ncol*sizeof(double));
  pomp_qr(y2,nrow,ncol,pivot,tau); // Q*R = Y*inv(diag(d))
  pomp_qr_x_inv_r(y2,nrow,ncol,y1,nrow); // y1 <- y1 %*% inv(R)

  // create Campbell weight vector
  d0 = sqrt(ncol)+alpha/sqrt(2.0);
  for (wbar = 0, i = 0; i < nrow; i++) {
    for (xx = 0, j = 0; j < ncol; j++) {
      x = y1[i+nrow*j];
      xx += x*x;
    }
    d = sqrt((nrow-1)*xx);	// Mahalonobis distance of row i
    if (d > d0) {
      x = d-d0;
      xx = exp(-0.5*x*x/beta)*d0/d;
    } else {
      xx = 1.0;
    }
    w[i] = xx;
    wbar += xx*xx;
  }
  wbar = sqrt(wbar-1);

  // compute weighted column means, center each column, precondition
  memcpy(y1,y,nrow*ncol*sizeof(double));
  for (yp = y1, j = 0; j < ncol; j++, yp += nrow) {
    for (xx = 0, x = 0, i = 0; i < nrow; i++) {
      xx += w[i]*yp[i];
      x += w[i];
    }
    xx /= x;			// column mean
    for (i = 0; i < nrow; i++) yp[i] -= xx; // center the column
    ydat[j] -= xx;		// subtract mean from realized probe
    for (xx = 0, i = 0; i < nrow; i++) {
      xx += yp[i]*yp[i];
      yp[i] /= wbar; 
    }
    x = sqrt(xx/(nrow-1)); // column SD
    for (i = 0; i < nrow; i++) yp[i] *= (w[i]/x); // precondition & weight
    ydat[j] /= x;
    half_log_det += log(x); // sum up logs(diag(D))
  }

  // do second QR decomposition & backsolve
  pomp_qr(y1,nrow,ncol,pivot,tau); // Q*R = diag(w)*Y*inv(diag(d))
  pomp_backsolve(y1,nrow,ncol,ydat,1,"Upper","Transpose","Non-unit"); // ydat <- ydat*inv(R)

  // compute residual sum of squares and add up logs of diag(R)
  xx = 0;
  for (yp = y1, rss = 0, i = nrow+1, j = 0; j < ncol; j++, yp += i) { // yp marches along the diagonal of R
    x = ydat[j];
    rss += x*x;
    half_log_det += log(fabs(*yp)); // log(diag(R))
    xx += log(fabs(*yp));
  }

  *loglik = -0.5*rss-half_log_det;
}


SEXP synth_loglik (SEXP ysim, SEXP ydat) {
  int nprotect = 0;
  SEXP loglik, dim, y;

  PROTECT(y = duplicate(AS_NUMERIC(ydat))); nprotect++;
  PROTECT(dim = GET_DIM(ysim)); nprotect++;
  PROTECT(ysim = AS_NUMERIC(ysim)); nprotect++;
  PROTECT(loglik = NEW_NUMERIC(1)); nprotect++;

  robust_synth_loglik(REAL(ysim),INTEGER(dim),REAL(y),REAL(loglik));

  UNPROTECT(nprotect);
  return loglik;
}

