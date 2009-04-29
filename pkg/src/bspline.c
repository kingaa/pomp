// -*- C++ -*-

#include "pomp_internal.h"

static void bspline_internal (double *y, const double *x, int nx, int i, int p, const double *knots, int nknots);

// B-spline basis

SEXP bspline_basis (SEXP x, SEXP degree, SEXP knots) {
  int nprotect = 0;
  SEXP y;
  int nx = length(x);
  int nknots = length(knots);
  int deg = INTEGER_VALUE(degree);
  int nbasis = nknots-deg-1;
  double *ydata;
  int i;
  if (deg < 0) error("must have degree > 0 in 'bspline.basis'");
  PROTECT(y = allocMatrix(REALSXP,nx,nbasis)); nprotect++;
  ydata = REAL(y);
  for (i = 0; i < nbasis; i++) {
    bspline_internal(ydata,REAL(x),nx,i,deg,REAL(knots),nknots);
    ydata += nx;
  }
  UNPROTECT(nprotect);
  return(y);
}

SEXP bspline_basis_function (SEXP x, SEXP i, SEXP degree, SEXP knots) {
  int nprotect = 0;
  SEXP y;
  int nx = length(x);
  int nknots = length(knots);
  int ival = INTEGER_VALUE(i);
  int deg = INTEGER_VALUE(degree);

  if ((ival < 0) || (deg < 0) || (ival + deg >= nknots-1)) 
    error("bad arguments in 'bspline'");

  PROTECT(y = NEW_NUMERIC(nx)); nprotect++;
  bspline_internal(REAL(y),REAL(x),nx,ival,deg,REAL(knots),nknots);
  UNPROTECT(nprotect);
  return(y);
}


static void bspline_internal (double *y, const double *x, int nx, int i, int p, const double *knots, int nknots)
{
  int j;
  double a, b;
  double y1[nx], y2[nx];
  int i2, p2;

  if (p == 0) {
    for (j = 0; j < nx; j++)
      y[j] = ((knots[i] <= x[j]) && (x[j] < knots[i+1]));
  } else {
    i2 = i+1;
    p2 = p-1;
    bspline_internal(y1,x,nx,i,p2,knots,nknots);
    bspline_internal(y2,x,nx,i2,p2,knots,nknots);
    for (j = 0; j < nx; j++) {
      a = (x[j]-knots[i]) / (knots[i+p]-knots[i]);
      b = (knots[i2+p]-x[j]) / (knots[i2+p]-knots[i2]);
      y[j] = a * y1[j] + b * y2[j];
    }
  }
}
