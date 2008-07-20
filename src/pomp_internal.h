// -*- C++ -*-

#ifndef _POMP_INTERNAL_H_
#define _POMP_INTERNAL_H_

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

#include "pomp.h"

# define MATCHROWNAMES(X,N) (matchnames(GET_ROWNAMES(GET_DIMNAMES(X)),(N)))
# define MATCHCOLNAMES(X,N) (matchnames(GET_COLNAMES(GET_DIMNAMES(X)),(N)))
# define MATCH_CHAR_TO_ROWNAMES(X,N,A) (match_char_to_names(GET_ROWNAMES(GET_DIMNAMES(X)),(N),(A)))

/* bspline.c */
SEXP bspline_basis(SEXP x, SEXP degree, SEXP knots);
SEXP bspline_basis_function(SEXP x, SEXP i, SEXP degree, SEXP knots);

/* dsobol.c */
SEXP sobol_sequence(SEXP dim);

/* pomp_fun.c */
SEXP pomp_fun_handler (SEXP pfun, int *use_native);

/* lookup_table.c */
SEXP lookup_in_table (SEXP ttable, SEXP xtable, SEXP t, int *index);

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

static SEXP matchnames (SEXP x, SEXP names) {
  int nprotect = 0;
  int n = length(names);
  int *idx, k;
  SEXP index, nm;
  PROTECT(nm = AS_CHARACTER(names)); nprotect++;
  PROTECT(index = match(x,names,0)); nprotect++;
  idx = INTEGER(index);
  for (k = 0; k < n; k++) {
    if (idx[k]==0) error("variable %s not found",CHARACTER_DATA(STRING_ELT(nm,k)));
    idx[k] -= 1;
  }
  UNPROTECT(nprotect);
  return index;
}

static SEXP match_char_to_names (SEXP x, int n, char **names) {
  int nprotect = 0;
  int *idx, k;
  SEXP index, nm;
  PROTECT(nm = NEW_CHARACTER(n)); nprotect++;
  for (k = 0; k < n; k++) {
    SET_STRING_ELT(nm,k,mkChar(names[k]));
  }
  PROTECT(index = match(x,nm,0)); nprotect++;
  idx = INTEGER(index);
  for (k = 0; k < n; k++) {
    if (idx[k]==0) {
      UNPROTECT(nprotect);
      error("variable %s not specified",names[k]);
    }
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
