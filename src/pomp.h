// -*- C++ -*-

#ifndef _POMP_H_
#define _POMP_H_

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>

#include "lookup_table.h"
#include "euler.h"

/* bspline.c */
SEXP bspline_basis(SEXP x, SEXP degree, SEXP knots);
SEXP bspline_basis_function(SEXP x, SEXP i, SEXP degree, SEXP knots);

/* dsobol.c */
SEXP sobol_sequence(SEXP dim);

/* resample.c */
SEXP systematic_resampling(SEXP weights);

static SEXP makearray (int rank, int *dim) {
  int nprotect = 0;
  int *dimp, k;
  SEXP dimx, x;
  PROTECT(dimx = NEW_INTEGER(rank)); nprotect++;
  dimp = INTEGER(dimx); 
  for (k = 0; k < rank; k++) dimp[k] = dim[k];
  PROTECT(x = allocArray(REALSXP,dimx)); nprotect++;
  UNPROTECT(nprotect);
  return x;
}

static SEXP matchrownames (SEXP x, SEXP names) {
  int nprotect = 0;
  int n = length(names);
  int *idx, k;
  SEXP index, nm;
  PROTECT(nm = AS_CHARACTER(names)); nprotect++;
  PROTECT(index = match(GET_ROWNAMES(GET_DIMNAMES(x)),names,0)); nprotect++;
  idx = INTEGER(index);
  for (k = 0; k < n; k++) {
    if (idx[k]==0) error("variable %s not specified",STRING_ELT(nm,k));
    idx[k] -= 1;
  }
  UNPROTECT(nprotect);
  return index;
}

static void setrownames (SEXP x, SEXP names, int n) {
  int nprotect = 0;
  SEXP dimnms, nm;
  PROTECT(nm = AS_CHARACTER(names)); nprotect++;
  PROTECT(dimnms = allocVector(VECSXP,n)); nprotect++;
  SET_ELEMENT(dimnms,0,nm);	// set row names
  SET_DIMNAMES(x,dimnms);
  UNPROTECT(nprotect);
}

static double expit (double x) {
  return 1.0/(1.0 + exp(-x));
}

static double logit (double x) {
  return log(x/(1-x));
}

#ifdef __cplusplus

template <class Scalar>
class view {
private:
  Scalar *data;
  int dim[2];
public:
  view (Scalar *x) {
    data = x;
    dim[0] = 0;
    dim[1] = 0;
  };
  view (Scalar *x, int d1) {
    data = x;
    dim[0] = d1;
    dim[1] = 0;
  };
  view (Scalar *x, int d1, int d2) {
    data = x;
    dim[0] = d1;
    dim[1] = d2;
  };
  ~view (void) {};
  inline Scalar& operator () (int d1) {
    return(data[d1]);
  };
  inline Scalar& operator () (int d1, int d2) {
    return(data[d1 + dim[0] * d2]);
  };
  inline Scalar& operator () (int d1, int d2, int d3) {
    return(data[d1 + dim[0] * (d2 + dim[1] * d3)]);
  };
};

#endif

#endif
