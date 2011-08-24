#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>

#include "pomp_internal.h"

SEXP lookup_in_table (SEXP ttable, SEXP xtable, SEXP t) {
  int nprotect = 0;
  int *dim, xdim[2], length, width;
  int j, nt;
  double *tp, *xp;
  SEXP X;

  PROTECT(t = AS_NUMERIC(t)); nprotect++;
  nt = LENGTH(t);

  dim = INTEGER(GET_DIM(xtable));
  length = dim[0]; width = dim[1];
  if (length != LENGTH(ttable)) {
    UNPROTECT(nprotect);
    error("incommensurate dimensions in 'lookup_in_table'");
  }

  if (nt > 1) {
    xdim[0] = width; xdim[1] = nt;
    PROTECT(X = makearray(2,xdim)); nprotect++;
    setrownames(X,GET_COLNAMES(GET_DIMNAMES(xtable)),2);
  } else {
    PROTECT(X = NEW_NUMERIC(width)); nprotect++;
    SET_NAMES(X,GET_COLNAMES(GET_DIMNAMES(xtable)));
  }

  struct lookup_table tab = {length,width,0,REAL(ttable),REAL(xtable)};
  for (j = 0, tp = REAL(t), xp = REAL(X); j < nt; j++, tp++, xp += width)
    table_lookup(&tab,*tp,xp,0);

  UNPROTECT(nprotect);
  return X;
}

// linear interpolation on a lookup table
void table_lookup (struct lookup_table *tab, double x, double *y, double *dydt)
{
  int flag = 0;
  int j, k;
  double e;
  tab->index = findInterval(tab->x,tab->length,x,TRUE,TRUE,tab->index,&flag);
  if (flag != 0)              // we are extrapolating
    warning("table_lookup: extrapolating (flag %d) at %le", flag, x);
  e = (x - tab->x[tab->index-1]) / (tab->x[tab->index] - tab->x[tab->index-1]);
  for (j = 0; j < tab->width; j++) {
    k = j*(tab->length)+(tab->index);
    y[j] = e*(tab->y[k])+(1-e)*(tab->y[k-1]);
    if (dydt != 0)
      dydt[j] = ((tab->y[k])-(tab->y[k-1]))/((tab->x[tab->index])-(tab->x[tab->index-1]));
  }
}

// compute the transmission coefficient using the basis functions
double dot_product (int dim, const double *basis, const double *coef)
{
  int j;
  double trans = 0.0;
  for (j = 0; j < dim; j++)
    trans += coef[j]*basis[j];
  return(trans);
}
