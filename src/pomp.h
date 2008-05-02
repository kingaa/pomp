// -*- C++ -*-

#ifndef _POMP_H_
#define _POMP_H_

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>

#include "interp.h"
#include "euler.h"

/* bspline.c */
SEXP bspline_basis(SEXP x, SEXP degree, SEXP knots);
SEXP bspline_basis_function(SEXP x, SEXP i, SEXP degree, SEXP knots);

/* dsobol.c */
SEXP sobol_sequence(SEXP dim);

/* resample.c */
SEXP systematic_resampling(SEXP weights);

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
