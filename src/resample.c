// -*- C++ -*-

#include <Rdefines.h>
#include "pomp_internal.h"

SEXP systematic_resampling(SEXP weights, SEXP np);
void nosort_resamp(int nw, double *w, int np, int *p, int offset);

SEXP systematic_resampling (SEXP weights, SEXP np)
{
  int m, n;
  SEXP perm;

  m = *(INTEGER(AS_INTEGER(np)));
  n = LENGTH(weights);
  PROTECT(perm = NEW_INTEGER(m));
  PROTECT(weights = AS_NUMERIC(weights));
  GetRNGstate();
  nosort_resamp(n,REAL(weights),m,INTEGER(perm),1);
  PutRNGstate();
  UNPROTECT(2);
  return(perm);
}

void nosort_resamp (int nw, double *w, int np, int *p, int offset)
{
  int i, j;
  double du, u;

  for (j = 1; j < nw; j++) w[j] += w[j-1];

  if (w[nw-1] <= 0.0)
    err("in 'systematic_resampling': non-positive sum of weights");

  du = w[nw-1] / ((double) np);
  u = -du*unif_rand();

  for (i = 0, j = 0; j < np; j++) {
    u += du;
    // In the following line, the second test is needed to correct
    // the infamous Bug of St. Patrick, 2017-03-17.
    while ((u > w[i]) && (i < nw-1)) i++;
    p[j] = i;
  }
  if (offset)			// add offset if needed
    for (j = 0; j < np; j++) p[j] += offset;

}
