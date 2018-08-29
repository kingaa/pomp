#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>

#include "pomp_internal.h"

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
  return tab;
}

SEXP lookup_in_table (SEXP covar, SEXP t) {
  int nprotect = 0;
  int xdim[2], nvar;
  int j, nt;
  double *tp, *xp;
  SEXP Cnames, X;

  PROTECT(t = AS_NUMERIC(t)); nprotect++;
  nt = LENGTH(t);
  PROTECT(Cnames = get_covariate_names(covar)); nprotect++;

  lookup_table_t tab = make_covariate_table(covar,&nvar);

  if (nt > 1) {
    xdim[0] = nvar; xdim[1] = nt;
    PROTECT(X = makearray(2,xdim)); nprotect++;
    setrownames(X,Cnames,2);
  } else {
    PROTECT(X = NEW_NUMERIC(nvar)); nprotect++;
    SET_NAMES(X,Cnames);
  }

  for (j = 0, tp = REAL(t), xp = REAL(X); j < nt; j++, tp++, xp += nvar)
    table_lookup(&tab,*tp,xp);

  UNPROTECT(nprotect);
  return X;
}

// linear interpolation on a lookup table
void table_lookup (lookup_table_t *tab, double x, double *y)
{
  int flag = 0;
  int j, k;
  double e;
  if ((tab == 0) || (tab->length < 1) || (tab->width < 1)) return;
  tab->index = findInterval(tab->x,tab->length,x,TRUE,TRUE,tab->index,&flag);
  // warn only if we are *outside* the interval
  if ((x < tab->x[0]) || (x > tab->x[(tab->length)-1]))
    warningcall(R_NilValue,"in 'table_lookup': extrapolating at %le.", x);
  e = (x - tab->x[tab->index-1]) / (tab->x[tab->index] - tab->x[tab->index-1]);
  for (j = 0; j < tab->width; j++) {
    k = j+(tab->width)*(tab->index);
    y[j] = e*(tab->y[k])+(1-e)*(tab->y[k-tab->width]);
    //    if (dydt != 0)
    //      dydt[j] = ((tab->y[k])-(tab->y[k-1]))/((tab->x[tab->index])-(tab->x[tab->index-1]));
  }
}

