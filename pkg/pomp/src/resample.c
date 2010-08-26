// -*- C++ -*-

#include "pomp_internal.h"
#include <Rdefines.h>

SEXP systematic_resampling (SEXP weights)
{
  int nprotect = 0;
  int n = length(weights);
  double u, du, *w;
  int i, j, *p;
  SEXP perm;

  PROTECT(perm = NEW_INTEGER(n)); nprotect++;

  p = INTEGER(perm);
  w = REAL(weights);
  for (j = 1; j < n; j++) w[j] += w[j-1];
  for (j = 0; j < n; j++) w[j] /= w[n-1];

  GetRNGstate();
  du = 1.0 / ((double) n);
  u = runif(-du,0);
  PutRNGstate();

  for (i = 0, j = 0; j < n; j++) {
    u += du;
    while (u > w[i]) i++;
    p[j] = i + 1; // must use 1-based indexing for compatibility with R!
  }

  UNPROTECT(nprotect);

  return(perm);
}

