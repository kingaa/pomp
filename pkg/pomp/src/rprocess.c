// dear emacs, please treat this as -*- C++ -*-

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <string.h>

#include "pomp_internal.h"

SEXP do_rprocess (SEXP object, SEXP xstart, SEXP times, SEXP params, SEXP offset, SEXP gnsi)
{
  int nprotect = 0;
  int *xdim, nvars, npars, nreps, nrepsx, ntimes, off;
  SEXP X, Xoff, copy, fn, fcall, rho;
  SEXP dimXstart, dimP, dimX;
  const char *dimnm[3] = {"variable","rep","time"};

  PROTECT(gnsi = duplicate(gnsi)); nprotect++;

  ntimes = length(times);
  if (ntimes < 2) {
    error("rprocess error: length(times)==0: no transitions, no work to do");
  }

  off = *(INTEGER(AS_INTEGER(offset)));
  if ((off < 0)||(off>=ntimes))
    error("illegal 'offset' value %d",off);

  PROTECT(xstart = as_matrix(xstart)); nprotect++;
  PROTECT(dimXstart = GET_DIM(xstart)); nprotect++;
  xdim = INTEGER(dimXstart);
  nvars = xdim[0]; nrepsx = xdim[1];

  PROTECT(params = as_matrix(params)); nprotect++;
  PROTECT(dimP = GET_DIM(params)); nprotect++;
  xdim = INTEGER(dimP);
  npars = xdim[0]; nreps = xdim[1]; 

  if (nrepsx > nreps) {		// more ICs than parameters
    if (nrepsx % nreps != 0) {
      error("rprocess error: larger number of replicates is not a multiple of smaller");
    } else {
      double *src, *tgt;
      int dims[2];
      int j, k;
      dims[0] = npars; dims[1] = nrepsx;
      PROTECT(copy = duplicate(params)); nprotect++;
      PROTECT(params = makearray(2,dims)); nprotect++;
      setrownames(params,GET_ROWNAMES(GET_DIMNAMES(copy)),2);
      src = REAL(copy);
      tgt = REAL(params);
      for (j = 0; j < nrepsx; j++) {
	for (k = 0; k < npars; k++, tgt++) {
	  *tgt = src[k+npars*(j%nreps)];
	}
      }
    }
    nreps = nrepsx;
  } else if (nrepsx < nreps) {	// more parameters than ICs
    if (nreps % nrepsx != 0) {
      error("rprocess error: larger number of replicates is not a multiple of smaller");
    } else {
      double *src, *tgt;
      int dims[2];
      int j, k;
      dims[0] = nvars; dims[1] = nreps;
      PROTECT(copy = duplicate(xstart)); nprotect++;
      PROTECT(xstart = makearray(2,dims)); nprotect++;
      setrownames(xstart,GET_ROWNAMES(GET_DIMNAMES(copy)),2);
      src = REAL(copy);
      tgt = REAL(xstart);
      for (j = 0; j < nreps; j++) {
	for (k = 0; k < nvars; k++, tgt++) {
	  *tgt = src[k+nvars*(j%nrepsx)];
	}
      }
    }
  }

  // extract the process function
  PROTECT(fn = GET_SLOT(object,install("rprocess"))); nprotect++;
  // construct the call
  PROTECT(fcall = VectorToPairList(GET_SLOT(object,install("userdata")))); nprotect++;
  PROTECT(fcall = LCONS(gnsi,fcall)); nprotect++;
  SET_TAG(fcall,install(".getnativesymbolinfo"));
  PROTECT(fcall = LCONS(GET_SLOT(object,install("zeronames")),fcall)); nprotect++;
  SET_TAG(fcall,install("zeronames"));
  PROTECT(fcall = LCONS(GET_SLOT(object,install("covar")),fcall)); nprotect++;
  SET_TAG(fcall,install("covar"));
  PROTECT(fcall = LCONS(GET_SLOT(object,install("tcovar")),fcall)); nprotect++;
  SET_TAG(fcall,install("tcovar"));
  PROTECT(fcall = LCONS(params,fcall)); nprotect++;
  SET_TAG(fcall,install("params"));
  PROTECT(fcall = LCONS(AS_NUMERIC(times),fcall)); nprotect++;
  SET_TAG(fcall,install("times"));
  PROTECT(fcall = LCONS(xstart,fcall)); nprotect++;
  SET_TAG(fcall,install("xstart"));
  PROTECT(fcall = LCONS(fn,fcall)); nprotect++;
  PROTECT(rho = (CLOENV(fn))); nprotect++; // environment of the function
  PROTECT(X = eval(fcall,rho)); nprotect++; // do the call
  PROTECT(dimX = GET_DIM(X)); nprotect++;
  if ((isNull(dimX)) || (length(dimX) != 3)) {
    error("rprocess error: user 'rprocess' must return a rank-3 array");
  }
  xdim = INTEGER(dimX);
  if ((xdim[0] != nvars) || (xdim[1] != nreps) || (xdim[2] != ntimes)) {
    error("rprocess error: user 'rprocess' must return a %d x %d x %d array",nvars,nreps,ntimes);
  }
  if (isNull(GET_ROWNAMES(GET_DIMNAMES(X)))) {
    error("rprocess error: user 'rprocess' must return an array with rownames");
  }
  if (off > 0) {
    xdim[2] -= off;
    PROTECT(Xoff = makearray(3,xdim)); nprotect++;
    setrownames(Xoff,GET_ROWNAMES(GET_DIMNAMES(X)),3);
    fixdimnames(Xoff,dimnm,3);
    memcpy(REAL(Xoff),REAL(X)+off*nvars*nreps,(ntimes-off)*nvars*nreps*sizeof(double));
    UNPROTECT(nprotect);
    return Xoff;
  } else {
    fixdimnames(X,dimnm,3);
    UNPROTECT(nprotect);
    return X;
  }
}
