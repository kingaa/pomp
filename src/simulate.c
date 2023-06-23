// -*- C++ -*-

#include <Rdefines.h>
#include "internal.h"

SEXP do_simulate (SEXP object, SEXP params, SEXP nsim, SEXP rettype, SEXP gnsi)
{

  SEXP t0, times, x0, x, y;
  SEXP ans = R_NilValue, ans_names = R_NilValue;
  SEXP simnames;
  int return_type = *(INTEGER(rettype)); // 0 = array, 1 = pomps

  if (LENGTH(nsim) != 1) err("'nsim' must be a single integer"); // #nocov

  PROTECT(params = as_matrix(params));

  PROTECT(t0 = GET_SLOT(object,install("t0")));
  PROTECT(times = GET_SLOT(object,install("times")));

  // initialize the simulations
  PROTECT(x0 = do_rinit(object,params,t0,nsim,gnsi));
  PROTECT(simnames = GET_COLNAMES(GET_DIMNAMES(x0)));

  // call 'rprocess' to simulate state process
  PROTECT(x = do_rprocess(object,x0,t0,times,params,gnsi));

  // call 'rmeasure' to simulate the measurement process
  PROTECT(y = do_rmeasure(object,x,times,params,gnsi));

  setcolnames(x,simnames);
  setcolnames(y,simnames);

  int nprotect = 7;

  switch (return_type) {

  case 0:  // return a list of arrays

    PROTECT(ans = NEW_LIST(2));
    PROTECT(ans_names = NEW_CHARACTER(2));
    nprotect += 2;
    SET_STRING_ELT(ans_names,0,mkChar("states"));
    SET_STRING_ELT(ans_names,1,mkChar("obs"));
    SET_NAMES(ans,ans_names);
    SET_ELEMENT(ans,0,x);
    SET_ELEMENT(ans,1,y);

    break;

  case 1: default:

    // a list to hold the pomp objects
    {
      SEXP pp, xx, yy, po;
      const int *xdim;
      int nvar, npar, nobs, nsim, ntim, nparsets;
      int dim[2], i, j, k;

      PROTECT(po = duplicate(object));
      SET_SLOT(po,install("t0"),t0);
      SET_SLOT(po,install("times"),times);

      xdim = INTEGER(GET_DIM(x));
      nvar = xdim[0]; nsim = xdim[1]; ntim = xdim[2];

      xdim = INTEGER(GET_DIM(y));
      nobs = xdim[0]; // second dimensions of 'x' and 'y' must agree

      xdim = INTEGER(GET_DIM(params));
      npar = xdim[0]; nparsets = xdim[1];

      dim[0] = nvar; dim[1] = ntim;
      PROTECT(xx = makearray(2,dim));
      setrownames(xx,GET_ROWNAMES(GET_DIMNAMES(x)),2);
      SET_SLOT(po,install("states"),xx);

      dim[0] = nobs; dim[1] = ntim;
      PROTECT(yy = makearray(2,dim));
      setrownames(yy,GET_ROWNAMES(GET_DIMNAMES(y)),2);
      SET_SLOT(po,install("data"),yy);

      PROTECT(pp = NEW_NUMERIC(npar));
      SET_NAMES(pp,GET_ROWNAMES(GET_DIMNAMES(params)));
      SET_SLOT(po,install("params"),pp);

      PROTECT(ans = NEW_LIST(nsim));
      SET_NAMES(ans,simnames);

      nprotect += 5;

      for (k = 0; k < nsim; k++) {

        SEXP po2;
        double *xs = REAL(x), *ys = REAL(y), *ps = REAL(params);
        double *xt, *yt, *pt;

        PROTECT(po2 = duplicate(po));
        xt = REAL(GET_SLOT(po2,install("states")));
        yt = REAL(GET_SLOT(po2,install("data")));
        pt = REAL(GET_SLOT(po2,install("params")));

        memcpy(pt,ps+npar*(k%nparsets),npar*sizeof(double));

        // copy x[,k,] and y[,k,] into po2
        for (j = 0; j < ntim; j++) {
          for (i = 0; i < nvar; i++, xt++) *xt = xs[i+nvar*(k+nsim*j)];
          for (i = 0; i < nobs; i++, yt++) *yt = ys[i+nobs*(k+nsim*j)];
        }

        SET_ELEMENT(ans,k,po2);
        UNPROTECT(1);

      }

      break;

    }

  }

  UNPROTECT(nprotect);
  return ans;

}
