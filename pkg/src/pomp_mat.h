#ifndef _POMP_MAT_H_
#define _POMP_MAT_H_

#include "pomp_internal.h"
#include <R_ext/Lapack.h>

static R_INLINE void pomp_backsolve (double *a, int m, int n, double *x, int incx) {
  // Level 2 BLAS triangular-matrix solver
  // DTRSV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
  F77_NAME(dtrsv)("U","N","N",&n,a,&m,x,&incx);
}

static R_INLINE void pomp_qr (double *a, int m, int n, int *pivot, double *tau) {
  int info, j, lwork = -1;
  double work1;
  for (j = 0; j < n; j++) pivot[j] = 0; // zero out the pivots (assumed by DGEQP3)
  // LAPACK QR decomposition routine
  // DGEQP3(M,N,A,LDA,JPVT,TAU,WORK,LWORK,INFO)
  F77_NAME(dgeqp3)(&m,&n,a,&m,pivot,tau,&work1,&lwork,&info); // workspace query
  lwork = (int) ceil(work1);
  {
    double work[lwork];
    F77_NAME(dgeqp3)(&m,&n,a,&m,pivot,tau,work,&lwork,&info); // actual call
  }
  for (j = 0; j < n; j++) pivot[j]--;
}

static R_INLINE void pomp_qrqy (double *c, double *a, double *tau, int m, int n, int k, int left, int tp) {
  int lda, info, lwork = -1;
  char side, trans;
  double work1;
 
  side = (left) ? 'L' : 'R';
  lda = (left) ? m : n;
  trans = (tp) ? 'T' : 'N';

  // workspace query
  // DORMQR(SIDE,TRANS,M,N,K,A,LDA,TAU,C,LDC,WORK,LWORK,INFO)
  F77_NAME(dormqr)(&side,&trans,&m,&n,&k,a,&lda,tau,c,&m,&work1,&lwork,&info);
  lwork = (int) ceil(work1);
  {				// actual call
    double work[lwork];
    F77_NAME(dormqr)(&side,&trans,&m,&n,&k,a,&lda,tau,c,&m,work,&lwork,&info);
  }
}

#endif
