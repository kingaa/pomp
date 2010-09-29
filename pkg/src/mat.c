// -*- mode: C++; -*-
// some functions for linear algebra written by Simon N. Wood.
// codes have been modified by AAK to suit his own peccadilloes

#include "pomp_internal.h"
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>

void pomp_backsolve (double *R, int *r, int *c, double *B, double *C, int *bc) {
  double x;
  int i, j, k;
  for (j = 0; j < *bc; j++) { // work across columns of B & C
    for (i = *c-1; i >= 0; i--) { // work up each column of B & C
      for (k = i+1, x = 0.0; k < *c; k++) 
	x += R[i+*r*k]*C[k+*c*j];
      C[i+*c*j] = (B[i+*c*j]-x)/R[i+*r*i];
    }
  }
}

void pomp_qr (double *x, int *r, int *c, int *pivot, double *tau) {
  int info, j, lwork = -1;
  double work1;
  // workspace query
  F77_NAME(dgeqp3)(r,c,x,r,pivot,tau,&work1,&lwork,&info);
  lwork = (int) floor(work1);
  if ((work1-lwork) > 0.5) lwork++;
  {
    double work[lwork];
    // actual call
    F77_NAME(dgeqp3)(r,c,x,r,pivot,tau,work,&lwork,&info); 
  }
  for (j = 0; j < *c; j++) pivot[j]--;
  // ... for 'tis C in which we work and not the 'cursed Fortran...
}


void pomp_qrqy (double *b, double *a, double *tau, int *r, int *c, int *k, int *left, int *tp) {
  char side, trans;
  int lda, info, lwork = -1;
  double work1;
 
  side = (*left) ? 'L' : 'R';
  lda = (*left) ? *r : *c;
  trans = (*tp) ? 'T' : 'N';
  // workspace query
  F77_NAME(dormqr)(&side,&trans,r,c,k,a,&lda,tau,b,r,&work1,&lwork,&info);
  lwork = (int) floor(work1);
  if ((work1-lwork) > 0.5) lwork++;
  {
    double work[lwork];
    // actual call
    F77_NAME(dormqr)(&side,&trans,r,c,k,a,&lda,tau,b,r,work,&lwork,&info); 
  }
}
