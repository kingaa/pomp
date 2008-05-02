#include "interp.h"

// linear interpolation on a lookup table
void table_lookup (struct lookup_table *tab, double x, double *y, double *dydt)
{
  int flag = 0;
  int j;
  double e;
  tab->index = findInterval(tab->x,tab->length,x,TRUE,TRUE,tab->index,&flag);
  if (flag != 0)              // we are extrapolating
    warning("table_lookup: extrapolating (flag %d) at %le", flag, x);
  e = (x - tab->x[tab->index-1]) / (tab->x[tab->index] - tab->x[tab->index-1]);
  for (j = 0; j < tab->width; j++) {
    y[j] = tab->y[j][tab->index] * e + tab->y[j][tab->index-1] * (1-e);
  }
  if (dydt != 0) {
    for (j = 0; j < tab->width; j++) {
      dydt[j] = (tab->y[j][tab->index] - tab->y[j][tab->index-1]) / (tab->x[tab->index] - tab->x[tab->index-1]);
    }
  }
}

// compute the transmission coefficient using the basis functions
// useful in combination with table_lookup
double dot_product (int dim, const double *basis, const double *coef)
{
  int j;
  double trans = 0.0;
  for (j = 0; j < dim; j++)
    trans += coef[j]*basis[j];
  return(trans);
}
