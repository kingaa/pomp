// -*- C++ -*-

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>

#include "pomp_internal.h"

SEXP LogitTransform(SEXP P);
SEXP ExpitTransform(SEXP X);
SEXP LogBarycentricTransform(SEXP X);
SEXP InverseLogBarycentricTransform(SEXP Y);

SEXP LogitTransform (SEXP P) {
  double *p;
  int k, n;
  PROTECT(P = duplicate(AS_NUMERIC(P)));
  p = REAL(P);
  n = LENGTH(P);
  for (k = 0; k < n; k++, p++)
    *p = logit(*p);
  UNPROTECT(1);
  return P;
}

SEXP ExpitTransform (SEXP X) {
  double *x;
  int k, n;
  PROTECT(X = duplicate(AS_NUMERIC(X)));
  x = REAL(X);
  n = LENGTH(X);
  for (k = 0; k < n; k++, x++)
    *x = 1.0/(1.0+exp(-(*x)));
  UNPROTECT(1);
  return X;
}

SEXP LogBarycentricTransform (SEXP X) {
  SEXP Y;
  PROTECT(X = AS_NUMERIC(X));
  PROTECT(Y = NEW_NUMERIC(LENGTH(X)));
  to_log_barycentric(REAL(Y),REAL(X),LENGTH(X));
  UNPROTECT(2);
  return Y;
}

SEXP InverseLogBarycentricTransform (SEXP Y) {
  SEXP X;
  PROTECT(Y = AS_NUMERIC(Y));
  PROTECT(X = NEW_NUMERIC(LENGTH(Y)));
  from_log_barycentric(REAL(X),REAL(Y),LENGTH(Y));
  UNPROTECT(2);
  return X;
}
