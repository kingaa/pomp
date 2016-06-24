#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>

#include "pomp_internal.h"

struct lookup_table make_covariate_table (SEXP object, int *ncovar) {
  struct lookup_table tab;
  int *dim;
  dim = INTEGER(GET_DIM(GET_SLOT(object,install("covar"))));
  tab.length = dim[0];
  *ncovar = tab.width = dim[1];
  tab.index = 0;
  tab.x = REAL(GET_SLOT(object,install("tcovar")));
  tab.y = REAL(GET_SLOT(object,install("covar")));
  return tab;
}

SEXP lookup_in_table (SEXP ttable, SEXP xtable, SEXP t) {
  int nprotect = 0;
  int *dim, xdim[2], ntimes, nvar;
  int j, nt;
  double *tp, *xp;
  SEXP X;

  PROTECT(t = AS_NUMERIC(t)); nprotect++;
  nt = LENGTH(t);

  dim = INTEGER(GET_DIM(xtable));
  ntimes = dim[0]; nvar = dim[1];
  PROTECT(ttable = AS_NUMERIC(ttable)); nprotect++;
  if (ntimes != LENGTH(ttable))
    errorcall(R_NilValue,"in 'lookup_in_table': incommensurate dimensions");

  if (nt > 1) {
    xdim[0] = nvar; xdim[1] = nt;
    PROTECT(X = makearray(2,xdim)); nprotect++;
    setrownames(X,GET_COLNAMES(GET_DIMNAMES(xtable)),2);
  } else {
    PROTECT(X = NEW_NUMERIC(nvar)); nprotect++;
    SET_NAMES(X,GET_COLNAMES(GET_DIMNAMES(xtable)));
  }

  struct lookup_table tab = {ntimes,nvar,0,REAL(ttable),REAL(xtable)};
  for (j = 0, tp = REAL(t), xp = REAL(X); j < nt; j++, tp++, xp += nvar)
    table_lookup(&tab,*tp,xp);

  UNPROTECT(nprotect);
  return X;
}

// linear interpolation on a lookup table
//void table_lookup (struct lookup_table *tab, double x, double *y, double *dydt)
void table_lookup (struct lookup_table *tab, double x, double *y)
{
  int flag = 0;
  int j, k;
  double e;
  if ((tab == 0) || (tab->length < 1) || (tab->width < 1)) return;
  tab->index = findInterval(tab->x,tab->length,x,TRUE,TRUE,tab->index,&flag);
  // warn only if we are *outside* the interval
  if ((x < tab->x[0]) || (x > tab->x[(tab->length)-1]))
    warningcall(R_NilValue,"in 'table_lookup': extrapolating at %le", x);
  e = (x - tab->x[tab->index-1]) / (tab->x[tab->index] - tab->x[tab->index-1]);
  for (j = 0; j < tab->width; j++) {
    k = j*(tab->length)+(tab->index);
    y[j] = e*(tab->y[k])+(1-e)*(tab->y[k-1]);
    //    if (dydt != 0)
    //      dydt[j] = ((tab->y[k])-(tab->y[k-1]))/((tab->x[tab->index])-(tab->x[tab->index-1]));
  }
}

