// -*- mode: C++; -*-
#include "internal.h"

SEXP apply_probe_data (SEXP object, SEXP probes) {
  SEXP retval, data, vals;
  int nprobe;
  int i;

  nprobe = LENGTH(probes);
  PROTECT(data = GET_SLOT(object,install("data")));
  PROTECT(vals = NEW_LIST(nprobe));
  SET_NAMES(vals,GET_NAMES(probes));

  for (i = 0; i < nprobe; i++) {
    SET_ELEMENT(vals,i,eval(PROTECT(lang2(VECTOR_ELT(probes,i),data)),
                            R_ClosureEnv(VECTOR_ELT(probes,i))));
    if (!IS_NUMERIC(VECTOR_ELT(vals,i))) {
      err("probe %d returns a non-numeric result",i+1);
    }
    UNPROTECT(1);
  }
  PROTECT(vals = VectorToPairList(vals));
  PROTECT(retval = eval(PROTECT(LCONS(install("c"),vals)),R_BaseEnv));

  UNPROTECT(5);
  return retval;
}

SEXP apply_probe_sim (SEXP object, SEXP nsim, SEXP params,
                      SEXP probes, SEXP datval, SEXP gnsi) {
  SEXP x, y, names, sims;
  SEXP returntype, retval, val, valnames;
  int nprobe, nsims, nobs, ntimes, nvals;
  int xdim[2];
  double *xp, *yp;
  int p, s, i, j, k, len0 = 0, len = 0;

  PROTECT(nsim = AS_INTEGER(nsim));
  if ((LENGTH(nsim)!=1) || (INTEGER(nsim)[0]<=0))
    err("'nsim' must be a positive integer."); // #nocov

  PROTECT(gnsi = duplicate(gnsi));

  // 'names' holds the names of the probe values
  // we get these from a previous call to 'apply_probe_data'
  nprobe = LENGTH(probes);
  nvals = LENGTH(datval);
  PROTECT(names = GET_NAMES(datval));
  PROTECT(returntype = NEW_INTEGER(1));
  *(INTEGER(returntype)) = 0;
  PROTECT(sims = do_simulate(object,params,nsim,returntype,gnsi));
  PROTECT(y = VECTOR_ELT(sims,1));
  *(INTEGER(gnsi)) = 0;

  nobs = INTEGER(GET_DIM(y))[0];
  nsims = INTEGER(GET_DIM(y))[1];
  ntimes = INTEGER(GET_DIM(y))[2];
  // set up temporary storage
  xdim[0] = nobs; xdim[1] = ntimes;
  PROTECT(x = makearray(2,xdim));
  setrownames(x,GET_ROWNAMES(GET_DIMNAMES(y)),2);

  // set up matrix to hold results
  xdim[0] = nsims; xdim[1] = nvals;
  PROTECT(retval = makearray(2,xdim));
  PROTECT(valnames = NEW_LIST(2));
  SET_ELEMENT(valnames,1,names);        // set column names
  SET_DIMNAMES(retval,valnames);

  for (p = 0, k = 0; p < nprobe; p++, k += len) { // loop over probes

    R_CheckUserInterrupt();

    for (s = 0; s < nsims; s++) { // loop over simulations

      // copy the data from y[,s,] to x[,]
      xp = REAL(x);
      yp = REAL(y)+nobs*s;
      for (j = 0; j < ntimes; j++, yp += nobs*nsims) {
        for (i = 0; i < nobs; i++, xp++) *xp = yp[i];
      }

      // evaluate the probe on the simulated data
      PROTECT(val = eval(PROTECT(lang2(VECTOR_ELT(probes,p),x)),
                         R_ClosureEnv(VECTOR_ELT(probes,p))));
      if (!IS_NUMERIC(val)) {
        err("probe %d returns a non-numeric result.",p+1);
      }

      len = LENGTH(val);
      if (s == 0)
        len0 = len;
      else if (len != len0) {
        err("variable-sized results returned by probe %d.",p+1);
      }
      if (k+len > nvals)
        err("probes return different number of values on different datasets.");

      xp = REAL(retval); yp = REAL(val);
      for (i = 0; i < len; i++) xp[s+nsims*(i+k)] = yp[i];

      UNPROTECT(2);
    }

  }
  if (k != nvals)
    err("probes return different number of values on different datasets."); // #nocov

  UNPROTECT(9);
  return retval;
}
