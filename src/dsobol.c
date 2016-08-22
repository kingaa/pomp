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
  double *dp;
  int flag[2], taus = 0;
  int k;
  F77_CALL(insobl)(flag, &d, &n, &taus); // error check and setup 
  if (!flag[0]) errorcall(R_NilValue,"dimension is too high");
  if (!flag[1]) errorcall(R_NilValue,"too many points requested");
  PROTECT(data = allocMatrix(REALSXP,d,n)); dp = REAL(data);
  for (k = 0; k < n; k++) F77_CALL(gosobl)(dp + k*d); // compute
  UNPROTECT(1);
  return(data);
}
