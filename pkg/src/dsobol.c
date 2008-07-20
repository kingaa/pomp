// -*- C++ -*-

#include "pomp_internal.h"

static void dsobol (double *data, int dim, int n);
void F77_NAME(insobl)(int *, int *, int *, int *);
void F77_NAME(gosobl)(double *);

SEXP sobol_sequence (SEXP dim)
{
  SEXP data;
  int d = INTEGER(dim)[0];
  int n = INTEGER(dim)[1];
  PROTECT(data = allocMatrix(REALSXP,d,n));
  dsobol(REAL(data),d,n);
  UNPROTECT(1);
  return(data);
}

static void dsobol (double *data, int dim, int n)
{
  int flag[2], taus = 0;
  int k;
  F77_CALL(insobl)(flag, &dim, &n, &taus);
  if (!flag[0]) error("dimension is not OK in 'dsobol'");
  if (!flag[1]) error("number of points requested is not OK in 'dsobol'");
  for (k = 0; k < n; k++) F77_CALL(gosobl)(data + k*dim);
}
