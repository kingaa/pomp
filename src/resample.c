// -*- C++ -*-

#include "pomp_internal.h"
#include <Rdefines.h>

void nosort_resamp (int nw, double *w, int np, int *p, int offset) 
{
  int i, j;
  double du, u;

  for (j = 1; j < nw; j++) w[j] += w[j-1];

  if (w[nw-1] <= 0.0)
    error("non-positive sum of weights");

  du = w[nw-1] / ((double) np);
  u = runif(-du,0);

  for (i = 0, j = 0; j < np; j++) {
    u += du;
    while (u > w[i]) i++;
    p[j] = i;
  }
  if (offset)			// add offset if needed
    for (j = 0; j < np; j++) p[j] += offset;

}

SEXP systematic_resampling (SEXP weights)
{
  int n;
  SEXP perm;

  n = LENGTH(weights);
  PROTECT(perm = NEW_INTEGER(n));
  GetRNGstate();
  nosort_resamp(n,REAL(weights),n,INTEGER(perm),1);
  PutRNGstate();
  UNPROTECT(1);
  return(perm);
}
