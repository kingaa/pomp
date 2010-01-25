// -*- C++ -*-

#include "pomp_internal.h"

void periodic_bspline_basis_eval (double x, double period, int degree, int nbasis, double *y);
static void bspline_internal (double *y, const double *x, int nx, int i, int p, const double *knots, int nknots);

// B-spline basis

SEXP bspline_basis (SEXP x, SEXP degree, SEXP knots) {
  int nprotect = 0;
  SEXP y, xr, kr;
  int nx = length(x);
  int nknots = length(knots);
  int deg = INTEGER_VALUE(degree);
  int nb;
  double *ydata;
  int i;
  nb = nknots-deg-1;
  if (deg < 0) error("bspline.basis error: must have degree > 0");
  if (nb<=deg) error("bspline.basis error: must have nbasis > degree");
  PROTECT(xr = AS_NUMERIC(x)); nprotect++;
  PROTECT(kr = AS_NUMERIC(knots)); nprotect++;
  PROTECT(y = allocMatrix(REALSXP,nx,nb)); nprotect++;
  ydata = REAL(y);
  for (i = 0; i < nb; i++) {
    bspline_internal(ydata,REAL(xr),nx,i,deg,REAL(kr),nknots);
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


SEXP periodic_bspline_basis (SEXP x, SEXP nbasis, SEXP degree, SEXP period) {
  int nprotect = 0;
  SEXP y, xr;
  int nx = length(x);
  int nb = INTEGER_VALUE(nbasis);
  int deg = INTEGER_VALUE(degree);
  double pd = NUMERIC_VALUE(period);
  int j, k;
  double *xrd, *ydata, *val;
  if (deg < 0) error("periodic_bspline_basis error: must have degree >= 0");
  if (nb <= 0) error("periodic_bspline_basis error: must have nbasis > 0");
  if (deg > nb) error("periodic_bspline_basis error: must have degree <= nbasis");
  if (pd <= 0.0) error("periodic_bspline_basis error: must have period > 0");
  PROTECT(xr = AS_NUMERIC(x)); nprotect++;
  xrd = REAL(xr);
  PROTECT(y = allocMatrix(REALSXP,nx,nb)); nprotect++;
  ydata = REAL(y);
  val = (double *) Calloc(nb,double);
  for (j = 0; j < nx; j++) {
    periodic_bspline_basis_eval(xrd[j],pd,deg,nb,val);
    for (k = 0; k < nb; k++) ydata[j+nx*k] = val[k];
  }
  Free(val);
  UNPROTECT(nprotect);
  return y;
}

void periodic_bspline_basis_eval (double x, double period, int degree, int nbasis, double *y)
{
  int nknots = nbasis+2*degree+1;
  int shift = (degree-1)/2;
  double knots[nknots];
  double yy[nknots];
  double dx;
  int j, k;
  if (period <= 0.0) error("periodic_bspline_basis_eval error: must have period > 0");
  if (nbasis <= 0) error("periodic_bspline_basis_eval error: must have nbasis > 0");
  if (degree < 0) error("periodic_bspline_basis_eval error: must have degree >= 0");
  if (nbasis < degree) error("periodic_bspline_basis_eval error: must have nbasis >= degree");
  dx = period/((double) nbasis);
  for (k = -degree; k <= nbasis+degree; k++) {
    knots[degree+k] = k*dx;
  }
  x = fmod(x,period);
  if (x < 0.0) x += period;
  for (k = 0; k < nknots; k++) {
    bspline_internal(&yy[k],&x,1,k,degree,&knots[0],nknots);
  }
  for (k = 0; k < degree; k++) yy[k] += yy[nbasis+k];
  for (k = 0; k < nbasis; k++) {
    j = (shift+k)%nbasis;
    y[k] = yy[j];
  }
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
