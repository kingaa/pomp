// -*- C++ -*-

#include "internal.h"

// Compute log(mean(exp(X))) accurately,
// optionally with one element dropped
SEXP logmeanexp (const SEXP X, const SEXP Drop) {
  int j, n = LENGTH(X);
  int k = *INTEGER(Drop)-1;     // zero-based index
  double *x = REAL(X);
  long double m = R_NegInf;
  long double s = 0;
  for (j = 0; j < n; j++) {
    if (j != k)
      m = (x[j] > m) ? (long double) x[j] : m;
  }
  for (j = 0; j < n; j++) {
    if (j != k)
      s += expl((long double) x[j] - m);
  }
  if (k >= 0 && k < n) n--;
  return ScalarReal(m + log(s/n));
}
