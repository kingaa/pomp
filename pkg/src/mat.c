// -*- mode: C++; -*-
// some functions for linear algebra written by Simon N. Wood.
// codes have been modified by AAK to suit his own peccadilloes

#include <math.h>
#include <R.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include "pomp_internal.h"

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
  int info, lwork = -1, i;
  double work1, *work;
  // workspace query
  F77_NAME(dgeqp3)(r,c,x,r,pivot,tau,&work1,&lwork,&info);
  lwork = (int) floor(work1);
  if ((work1-lwork) >0.5) lwork++;
  work = (double *) Calloc(lwork,double);
  // actual call
  F77_NAME(dgeqp3)(r,c,x,r,pivot,tau,work,&lwork,&info); 
  Free(work);
  for (i = 0; i < *c; i++) pivot[i]--;
  // ... for 'tis C in which we work and not the 'cursed Fortran...
}


void pomp_qrqy (double *b, double *a, double *tau, int *r, int *c, int *k, int *left, int *tp) {
  char side, trans='N';
  int lda, lwork = -1, info;
  double *work, work1;
 
  if (! *left) { 
    side='R';
    lda = *c;
  } else {
    side = 'L';
    lda= *r;
  }
  if (*tp) trans='T'; 
  // workspace query
  F77_NAME(dormqr)(&side,&trans,r,c,k,a,&lda,tau,b,r,&work1,&lwork,&info);
  lwork = (int) floor(work1);
  if ((work1-lwork) > 0.5) lwork++;
  work = (double *) Calloc(lwork,double);
  // actual call
  F77_NAME(dormqr)(&side,&trans,r,c,k,a,&lda,tau,b,r,work,&lwork,&info); 
  Free(work);
}
