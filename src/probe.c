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

SEXP apply_probe_sim (SEXP object, SEXP nsim, SEXP params, SEXP seed, SEXP probes, SEXP datval) {
  int nprotect = 0;
  SEXP y, obs, times, t0, call, names;
  SEXP retval, val, valnames, x;
  int nprobe, nsims, nvars, ntimes, nvals;
  int xdim[2];
  double *xp, *yp;
  int p, s, i, j, k, len0 = 0, len = 0;

  PROTECT(nsim = AS_INTEGER(nsim)); nprotect++;
  if ((LENGTH(nsim)!=1) || (INTEGER(nsim)[0]<=0))
    errorcall(R_NilValue,"'nsim' must be a positive integer");

  // 'names' holds the names of the probe values
  // we get these from a previous call to 'apply_probe_data'
  nprobe = LENGTH(probes);
  nvals = LENGTH(datval);
  PROTECT(names = GET_NAMES(datval)); nprotect++;
  PROTECT(t0 = GET_SLOT(object,install("t0"))); nprotect++;
  PROTECT(times = GET_SLOT(object,install("times"))); nprotect++;

  // call 'simulate' to get simulated data sets
  PROTECT(obs = NEW_LOGICAL(1)); nprotect++;
  LOGICAL(obs)[0] = 1;		// we set obs=TRUE
  PROTECT(call = LCONS(t0,R_NilValue)); nprotect++;
  SET_TAG(call,install("t0"));
  PROTECT(call = LCONS(times,call)); nprotect++;
  SET_TAG(call,install("times"));
  PROTECT(call = LCONS(obs,call)); nprotect++;
  SET_TAG(call,install("obs"));
  PROTECT(call = LCONS(params,call)); nprotect++;
  SET_TAG(call,install("params"));
  PROTECT(call = LCONS(seed,call)); nprotect++;
  SET_TAG(call,install("seed"));
  PROTECT(call = LCONS(nsim,call)); nprotect++;
  SET_TAG(call,install("nsim"));
  PROTECT(call = LCONS(object,call)); nprotect++;
  SET_TAG(call,install("object"));
  PROTECT(call = LCONS(install("simulate"),call)); nprotect++;
  PROTECT(y = eval(call,R_GlobalEnv)); nprotect++;

  nvars = INTEGER(GET_DIM(y))[0];
  nsims = INTEGER(GET_DIM(y))[1];
  ntimes = INTEGER(GET_DIM(y))[2]; // recall that 'simulate' returns a value for time zero

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
