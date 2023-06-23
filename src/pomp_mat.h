#ifndef _POMP_MAT_H_
#define _POMP_MAT_H_

#include <Rconfig.h>
#include <R_ext/Lapack.h>
#ifndef FCONE
# define FCONE
#endif
#include "internal.h"

static R_INLINE void pomp_backsolve (double *a, int lda, int n, double *x, int incx,
				     char *uplo, char *transpose, char *unit) {
  // Level 2 BLAS triangular-matrix solver DTRSV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
  // DTRSV:  x <- A ^{-1} x -or- x <- A ^{-T} x, A triangular
  // N is the order of A
  // LDA is A's leading dimension
  // INCX is the increment between successive X locations
  // UPLO is "U" or "L" depending on whether A is upper or lower triangular
  // TRANSPOSE is "T" or "N" depending on whether the transpose is desired
  // DIAG is "U" or "N" depending on whether A is unit triangular or not
  F77_CALL(dtrsv)(uplo,transpose,unit,&n,a,&lda,x,&incx FCONE FCONE FCONE);
}

static R_INLINE void pomp_qr (double *a, int m, int n, int *pivot, double *tau) {
  int info, j, lwork = -1;
  double work1;
  for (j = 0; j < n; j++) pivot[j] = 0; // zero out the pivots (assumed by DGEQP3)
  // LAPACK QR decomposition routine
  // DGEQP3(M,N,A,LDA,JPVT,TAU,WORK,LWORK,INFO)
  F77_CALL(dgeqp3)(&m,&n,a,&m,pivot,tau,&work1,&lwork,&info); // workspace query
  lwork = (int) ceil(work1);
  {
    double work[lwork];
    F77_CALL(dgeqp3)(&m,&n,a,&m,pivot,tau,work,&lwork,&info); // actual call
  }
  for (j = 0; j < n; j++) pivot[j]--;
}

static R_INLINE void pomp_qrqy (double *c, double *a, int lda, double *tau, int m, int n, int k, char *side, char *transpose) {
  int info, lwork = -1;
  double work1;
 
  // workspace query
  // DORMQR(SIDE,TRANS,M,N,K,A,LDA,TAU,C,LDC,WORK,LWORK,INFO)
  F77_CALL(dormqr)(side,transpose,&m,&n,&k,a,&lda,tau,c,&m,&work1,&lwork,&info FCONE FCONE);
  lwork = (int) ceil(work1);
  {				// actual call
    double work[lwork];
    F77_CALL(dormqr)(side,transpose,&m,&n,&k,a,&lda,tau,c,&m,work,&lwork,&info FCONE FCONE);
  }
}

#endif
