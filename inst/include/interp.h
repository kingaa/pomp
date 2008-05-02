// -*- C++ -*-

#ifndef _POMP_INTERP_H_
#define _POMP_INTERP_H_

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>

struct lookup_table {
  int length, width;
  int index;
  double *x;
  double **y;
};

double dot_product (int dim, const double *basis, const double *coef);
void table_lookup (struct lookup_table *tab, double x, double *y, double *dydt);

#endif
