// -*- C++ -*-

#include "pomp_internal.h"

void periodic_bspline_basis_eval (double x, double period, int degree, int nbasis, double *y);
static void bspline_internal (double *y, const double *x, int nx, int i, int p, const double *knots);

// B-spline basis

SEXP bspline_basis (SEXP x, SEXP nbasis, SEXP degree) {
  int nprotect = 0;
  SEXP y, xr;
  int nx = length(x);
  int nb = INTEGER_VALUE(nbasis);
  int deg = INTEGER_VALUE(degree);
  int nk = nb+deg+1;
  double dx, minx, maxx;
  double knots[nk];
  double *xdata, *ydata;
  int i;
  if (deg < 0) error("bspline.basis error: must have degree > 0");
  if (nb <= deg) error("bspline.basis error: must have nbasis > degree");
  PROTECT(xr = AS_NUMERIC(x)); nprotect++;
  PROTECT(y = allocMatrix(REALSXP,nx,nb)); nprotect++;
  xdata = REAL(xr);
  ydata = REAL(y);
  for (i = 1, minx = maxx = xdata[0]; i < nx; i++) {
    minx = (minx > xdata[i]) ? xdata[i] : minx;
    maxx = (maxx < xdata[i]) ? xdata[i] : maxx;
  }
  dx = (maxx-minx)/((double) (nb-deg));
  knots[0] = minx-deg*dx;
  for (i = 1; i < nk; i++) knots[i] = knots[i-1]+dx;
  for (i = 0; i < nb; i++) {
    bspline_internal(ydata,xdata,nx,i,deg,&knots[0]);
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
  bspline_internal(REAL(y),REAL(x),nx,ival,deg,REAL(knots));
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
  if (nb < deg) error("periodic_bspline_basis error: must have nbasis >= degree");
  if (pd <= 0.0) error("periodic_bspline_basis error: must have period > 0");
  PROTECT(xr = AS_NUMERIC(x)); nprotect++;
  xrd = REAL(xr);
  PROTECT(y = allocMatrix(REALSXP,nx,nb)); nprotect++;
  ydata = REAL(y);
  val = (double *) R_alloc(nb,sizeof(double));
  for (j = 0; j < nx; j++) {
    periodic_bspline_basis_eval(xrd[j],pd,deg,nb,val);
    for (k = 0; k < nb; k++) ydata[j+nx*k] = val[k];
  }
  UNPROTECT(nprotect);
  return y;
}

void periodic_bspline_basis_eval (double x, double period, int degree, int nbasis, double *y)
{
  int nknots = nbasis+2*degree+1;
  int shift = (degree-1)/2;
  double *knots = NULL, *yy = NULL;
  double dx;
  int j, k;

  if (period <= 0.0) error("periodic_bspline_basis_eval error: must have period > 0");
  if (nbasis <= 0) error("periodic_bspline_basis_eval error: must have nbasis > 0");
  if (degree < 0) error("periodic_bspline_basis_eval error: must have degree >= 0");
  if (nbasis < degree) error("periodic_bspline_basis_eval error: must have nbasis >= degree");

  knots = (double *) Calloc(nknots+degree+1,double);
  yy = (double *) Calloc(nknots,double);

  dx = period/((double) nbasis);
  for (k = -degree; k <= nbasis+degree; k++) {
    knots[degree+k] = k*dx;
  }
  x = fmod(x,period);
  if (x < 0.0) x += period;
  for (k = 0; k < nknots; k++) {
    bspline_internal(&yy[k],&x,1,k,degree,knots);
  }
  for (k = 0; k < degree; k++) yy[k] += yy[nbasis+k];
  for (k = 0; k < nbasis; k++) {
    j = (shift+k)%nbasis;
    y[k] = yy[j];
  }
  Free(yy); Free(knots);
}

static void bspline_internal (double *y, const double *x, int nx, int i, int p, const double *knots)
{
  int j;
  if (p == 0) {
    for (j = 0; j < nx; j++)
      y[j] = (double) ((knots[i] <= x[j]) && (x[j] < knots[i+1]));
  } else {
    int i2 = i+1;
    int p2 = p-1;
    double *y1 = (double *) Calloc(nx,double);
    double *y2 = (double *) Calloc(nx,double);
    double a, b;
    bspline_internal(y1,x,nx,i,p2,knots);
    bspline_internal(y2,x,nx,i2,p2,knots);
    for (j = 0; j < nx; j++) {
      a = (x[j]-knots[i]) / (knots[i+p]-knots[i]);
      b = (knots[i2+p]-x[j]) / (knots[i2+p]-knots[i2]);
      y[j] = a * y1[j] + b * y2[j];
    }
    Free(y1); Free(y2);
  }
}
