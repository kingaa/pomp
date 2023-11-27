// -*- mode: C++ -*-

#include "internal.h"
#include "pomp_mat.h"

// Campbell's robust variance-covariance estimation approach
// described on p. 231 of Krzanowski 1988
// with additional pre-conditioning for numerical stability
// translated into C from R code by Simon N. Wood

static void robust_synth_loglik (double *y, int *dim, double *ydat, double *loglik) {
  // 'ydat' is destroyed
  // 'y' is preserved
  int nrow = dim[0];
  int ncol = dim[1];
  double w[nrow], tau[ncol], work[ncol];
  int info = 0;
  double one = 1.0;
  double *y1, *y2, *yp;
  long double x, xx, wbar, d, d0, rss, half_log_det;
  long double alpha = 2.0, beta = 1.25;
  int i, j;

  half_log_det = ncol*M_LN_SQRT_2PI;

  y1 = (double *) R_alloc(nrow*ncol,sizeof(double));
  y2 = (double *) R_alloc(nrow*ncol,sizeof(double));

  if (nrow <= ncol)
    err("'nsim' (=%d) should be (much) larger than the number of probes (=%d)",nrow,ncol); // #nocov

  // compute column means, center each column, precondition
  memcpy(y1,y,nrow*ncol*sizeof(double));
  for (yp = y1, j = 0; j < ncol; j++, yp += nrow) {
    for (x = 0, i = 0; i < nrow; i++) x += yp[i];
    x /= nrow;
    for (i = 0; i < nrow; i++) yp[i] -= x; // center the column
    for (x = 0, i = 0; i < nrow; i++) x += yp[i]*yp[i];
    d = sqrt(x/(nrow-1));		   // column SD
    for (i = 0; i < nrow; i++) yp[i] /= d; // precondition
  }

  // do first QR decomposition & backsolve
  memcpy(y2,y1,nrow*ncol*sizeof(double));
  // LAPACK QR decomposition without pivoting DGEQR2(M,N,A,LDA,TAU,WORK,INFO)
  F77_CALL(dgeqr2)(&nrow,&ncol,y2,&nrow,tau,work,&info);
  // Level-3 BLAS triangular matrix solver DTRSM(SIDE,UPLO,TRANS,DIAG,M,N,ALPHA,A,LDA,B,LDB)
  F77_CALL(dtrsm)("right","upper","no transpose","non-unit",&nrow,&ncol,&one,y2,&nrow,y1,&nrow FCONE FCONE FCONE FCONE);

  // create Campbell weight vector
  d0 = sqrt(ncol)+alpha/sqrt(2.0);
  for (wbar = 0, i = 0; i < nrow; i++) {
    for (xx = 0, j = 0; j < ncol; j++) {
      x = y1[i+nrow*j];
      xx += x*x;
    }
    d = sqrt((nrow-1)*xx); // Mahalonobis distance of row i
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
    for (x = 0, xx = 0, i = 0; i < nrow; i++) {
      x += w[i];
      xx += w[i]*yp[i];
    }
    xx /= x;			// column mean
    for (i = 0; i < nrow; i++) yp[i] -= xx; // center the column
    ydat[j] -= xx;		// subtract mean from realized probe
    for (xx = 0, i = 0; i < nrow; i++) {
      xx += yp[i]*yp[i];
      yp[i] /= wbar;
    }
    d = sqrt(xx/(nrow-1)); // column SD
    for (i = 0; i < nrow; i++) yp[i] *= (w[i]/d); // precondition & weight
    ydat[j] /= d;
    half_log_det += log(d); // sum up logs(diag(D))
  }

  // do second QR decomposition & backsolve
  // LAPACK QR decomposition without pivoting DGEQR2(M,N,A,LDA,TAU,WORK,INFO)
  F77_CALL(dgeqr2)(&nrow,&ncol,y1,&nrow,tau,work,&info);
  pomp_backsolve(y1,nrow,ncol,ydat,1,"upper","transpose","non-unit");

  // compute residual sum of squares and add up logs of diag(R)
  for (yp = y1, rss = 0, i = nrow+1, j = 0; j < ncol; j++, yp += i) { // yp marches along the diagonal of R
    x = ydat[j];
    rss += x*x;
    half_log_det += log(fabs(*yp)); // log(diag(R))
  }

  *loglik = -0.5*rss-half_log_det;
}

SEXP synth_loglik (SEXP ysim, SEXP ydat) {
  SEXP loglik, dim, y;

  PROTECT(y = duplicate(AS_NUMERIC(ydat)));
  PROTECT(dim = GET_DIM(ysim));
  PROTECT(ysim = AS_NUMERIC(ysim));
  PROTECT(loglik = NEW_NUMERIC(1));

  robust_synth_loglik(REAL(ysim),INTEGER(dim),REAL(y),REAL(loglik));

  UNPROTECT(4);
  return loglik;
}
