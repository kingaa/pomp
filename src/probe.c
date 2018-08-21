// -*- mode: C++; -*-
#include "pomp_internal.h"

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
      CLOENV(VECTOR_ELT(probes,i))));
    if (!IS_NUMERIC(VECTOR_ELT(vals,i))) {
      errorcall(R_NilValue,"probe %ld returns a non-numeric result",i+1);
    }
    UNPROTECT(1);
  }
  PROTECT(vals = VectorToPairList(vals));
  PROTECT(retval = eval(PROTECT(LCONS(install("c"),vals)),R_BaseEnv));

  UNPROTECT(5);
  return retval;
}

SEXP apply_probe_sim (SEXP object, SEXP nsim, SEXP params, SEXP seed,
  SEXP probes, SEXP datval, SEXP gnsi) {
  int nprotect = 0;
  SEXP x0, x, y, times, t0, timesplus, offset, names;
  SEXP retval, val, valnames;
  int nprobe, nsims, nvars, ntimes, nvals;
  int xdim[2];
  double *xp, *yp;
  int p, s, i, j, k, len0 = 0, len = 0;

  PROTECT(nsim = AS_INTEGER(nsim)); nprotect++;
  if ((LENGTH(nsim)!=1) || (INTEGER(nsim)[0]<=0))
    errorcall(R_NilValue,"'nsim' must be a positive integer");

  PROTECT(gnsi = duplicate(gnsi)); nprotect++;

  // 'names' holds the names of the probe values
  // we get these from a previous call to 'apply_probe_data'
  nprobe = LENGTH(probes);
  nvals = LENGTH(datval);
  PROTECT(names = GET_NAMES(datval)); nprotect++;
  PROTECT(t0 = GET_SLOT(object,install("t0"))); nprotect++;
  PROTECT(times = GET_SLOT(object,install("times"))); nprotect++;
  // generate simulated data sets
  // first, the initial state
  PROTECT(x0 = do_rinit(object,params,t0,nsim,gnsi)); nprotect++;
  // now the 'rprocess'
  // recall that 'rprocess' needs to have 't0' in the 'times' vector
  // and we compensate for this with 'offset'
  PROTECT(timesplus = eval(lang3(install("c"),t0,times),R_BaseEnv)); nprotect++;
  PROTECT(offset = NEW_INTEGER(1)); nprotect++;
  *(INTEGER(offset)) = 1;
  PROTECT(x = do_rprocess(object,x0,timesplus,params,offset,gnsi)); nprotect++;
  // now the 'rmeasure'
  PROTECT(y = do_rmeasure(object,x,times,params,gnsi)); nprotect++;
  *(INTEGER(gnsi)) = 0;

  nvars = INTEGER(GET_DIM(y))[0];
  nsims = INTEGER(GET_DIM(y))[1];
  ntimes = INTEGER(GET_DIM(y))[2];
  // set up temporary storage
  xdim[0] = nvars; xdim[1] = ntimes;
  PROTECT(x = makearray(2,xdim)); nprotect++;
  setrownames(x,GET_ROWNAMES(GET_DIMNAMES(y)),2);

  // set up matrix to hold results
  xdim[0] = nsims; xdim[1] = nvals;
  PROTECT(retval = makearray(2,xdim)); nprotect++;
  PROTECT(valnames = NEW_LIST(2)); nprotect++;
  SET_ELEMENT(valnames,1,names);	// set column names
  SET_DIMNAMES(retval,valnames);

  for (p = 0, k = 0; p < nprobe; p++, k += len) { // loop over probes

    R_CheckUserInterrupt();

    for (s = 0; s < nsims; s++) { // loop over simulations

      // copy the data from y[,s,] to x[,]
      xp = REAL(x);
      yp = REAL(y)+nvars*s;
      for (j = 0; j < ntimes; j++, yp += nvars*nsims) {
        for (i = 0; i < nvars; i++, xp++) *xp = yp[i];
      }

      // evaluate the probe on the simulated data
      PROTECT(val = eval(PROTECT(lang2(VECTOR_ELT(probes,p),x)),
        CLOENV(VECTOR_ELT(probes,p))));
      if (!IS_NUMERIC(val)) {
        errorcall(R_NilValue,"probe %ld returns a non-numeric result",p+1);
      }

      len = LENGTH(val);
      if (s == 0)
        len0 = len;
      else if (len != len0) {
        errorcall(R_NilValue,"variable-sized results returned by probe %ld",p+1);
      }
      if (k+len > nvals)
        errorcall(R_NilValue,"probes return different number of values on different datasets");

      xp = REAL(retval); yp = REAL(val);
      for (i = 0; i < len; i++) xp[s+nsims*(i+k)] = yp[i];

      UNPROTECT(2);
    }

  }
  if (k != nvals)
    errorcall(R_NilValue,"probes return different number of values on different datasets");

  UNPROTECT(nprotect);
  return retval;
}
