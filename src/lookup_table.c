#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>

#include "internal.h"

SEXP get_covariate_names (SEXP object) {
  return GET_ROWNAMES(GET_DIMNAMES(GET_SLOT(object,install("table"))));
}

lookup_table_t make_covariate_table (SEXP object, int *ncovar) {
  lookup_table_t tab;
  int *dim;
  dim = INTEGER(GET_DIM(GET_SLOT(object,install("table"))));
  *ncovar = tab.width = dim[0];
  tab.length = dim[1];
  tab.index = 0;
  tab.x = REAL(GET_SLOT(object,install("times")));
  tab.y = REAL(GET_SLOT(object,install("table")));
  tab.order = *(INTEGER(GET_SLOT(object,install("order"))));
  return tab;
}

SEXP lookup_in_table (SEXP covar, SEXP t) {
  int xdim[2], nvar;
  int j, nt;
  double *tp, *xp;
  SEXP Cnames, X;

  PROTECT(t = AS_NUMERIC(t));
  nt = LENGTH(t);
  PROTECT(Cnames = get_covariate_names(covar));

  lookup_table_t tab = make_covariate_table(covar,&nvar);

  if (nt > 1) {
    xdim[0] = nvar; xdim[1] = nt;
    PROTECT(X = makearray(2,xdim));
    setrownames(X,Cnames,2);
  } else {
    PROTECT(X = NEW_NUMERIC(nvar));
    SET_NAMES(X,Cnames);
  }

  for (j = 0, tp = REAL(t), xp = REAL(X); j < nt; j++, tp++, xp += nvar)
    table_lookup(&tab,*tp,xp);

  UNPROTECT(3);
  return X;
}

// linear interpolation on a lookup table
void table_lookup (lookup_table_t *tab, double x, double *y)
{
  int flag = 0;
  int j, k, n;
  double e;
  if ((tab == 0) || (tab->length < 1) || (tab->width < 1)) return;
  tab->index = findInterval(tab->x,tab->length,x,1,1,tab->index,&flag);
  // warn only if we are *outside* the interval
  if ((x < tab->x[0]) || (x > tab->x[tab->length-1]))
    warn("in 'table_lookup': extrapolating at %le.", x);
  switch (tab->order) {
  case 1: default: // linear interpolation
    e = (x - tab->x[tab->index-1]) / (tab->x[tab->index] - tab->x[tab->index-1]);
    for (j = 0, k = tab->index*tab->width, n = k-tab->width; j < tab->width; j++, k++, n++) {
      y[j] = e*(tab->y[k])+(1-e)*(tab->y[n]);
    }
    break;
  case 0: // piecewise constant
    if (flag < 0) n = 0;
    else if (flag > 0) n = tab->index;
    else n = tab->index-1;
    for (j = 0, k = n*tab->width; j < tab->width; j++, k++) {
      y[j] = tab->y[k];
    }
    break;
  }
}
