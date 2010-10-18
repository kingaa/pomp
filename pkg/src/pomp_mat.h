#ifndef _POMP_MAT_H_
#define _POMP_MAT_H_

#include "pomp_internal.h"
#include <R_ext/Lapack.h>

static R_INLINE void pomp_backsolve (double *a, int m, int n, double *x, int incx, 
				     char *uplo, char *transpose, char *unit) {
  // Level 2 BLAS triangular-matrix solver
  // DTRSV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
  F77_NAME(dtrsv)(uplo,transpose,unit,&n,a,&m,x,&incx);
}

static R_INLINE void pomp_qr_x_inv_r (double *a, int m, int n, double *b, int ldb) {
  // Level 3 BLAS triangular-matrix solver
  // DTRSM(SIDE,UPLO,TRANS,DIAG,M,N,ALPHA,A,LDA,B,LDB)
  double alpha = 1.0;
  F77_NAME(dtrsm)("Right","Upper","No transpose","Non-unit",&m,&n,&alpha,a,&m,b,&ldb);
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
  char *side, *trans;
  double work1;
 
  side = (left) ? "left" : "right";
  lda = (left) ? m : n;
  trans = (tp) ? "transpose" : "no transpose";

  // workspace query
  // DORMQR(SIDE,TRANS,M,N,K,A,LDA,TAU,C,LDC,WORK,LWORK,INFO)
  F77_NAME(dormqr)(side,trans,&m,&n,&k,a,&lda,tau,c,&m,&work1,&lwork,&info);
  lwork = (int) ceil(work1);
  {				// actual call
    double work[lwork];
    F77_NAME(dormqr)(side,trans,&m,&n,&k,a,&lda,tau,c,&m,work,&lwork,&info);
  }
}

#endif
