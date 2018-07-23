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
  int *xdim, type, nvars, npars, nreps, nrepsx, ntimes, off;
  SEXP X, Xoff, copy, rproc, args, zeronames, tcovar, covar;
  SEXP dimXstart, dimP, dimX;
  const char *dimnm[3] = {"variable","rep","time"};

  PROTECT(gnsi = duplicate(gnsi)); nprotect++;

  PROTECT(times = AS_NUMERIC(times)); nprotect++;
  ntimes = length(times);
  if (ntimes < 2) {
    errorcall(R_NilValue,"in 'rprocess': length(times) < 2: with no transitions, there is no work to do.");
  }

  off = *(INTEGER(AS_INTEGER(offset)));
  if ((off < 0)||(off>=ntimes))
    errorcall(R_NilValue,"illegal 'offset' value %d",off);

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
      errorcall(R_NilValue,"in 'rprocess': the larger number of replicates is not a multiple of smaller.");
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
      errorcall(R_NilValue,"in 'rprocess': the larger number of replicates is not a multiple of smaller.");
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

  PROTECT(rproc = GET_SLOT(object,install("rprocess"))); nprotect++;
  PROTECT(args = VectorToPairList(GET_SLOT(object,install("userdata")))); nprotect++;
  PROTECT(zeronames = GET_SLOT(object,install("zeronames"))); nprotect++;
  PROTECT(covar = GET_SLOT(object,install("covar"))); nprotect++;
  PROTECT(tcovar = GET_SLOT(object,install("tcovar"))); nprotect++;

  // extract the process function
  type = *(INTEGER(GET_SLOT(rproc,install("type"))));
  switch (type) {
  case 2: // one-step simulator
  {
    SEXP fn;
    double deltat = 1.0;
    int method = 1;
    PROTECT(fn = GET_SLOT(rproc,install("step.fn"))); nprotect++;
    PROTECT(X = euler_model_simulator(
      fn,xstart,times,params,deltat,method,
      zeronames,tcovar,covar,args,gnsi)); nprotect++;
  }
    break;
  case 3: case 4: // discrete-time and Euler
  {
    SEXP fn;
    double deltat;
    int method = (type == 3) ? 2 : 0;
    PROTECT(fn = GET_SLOT(rproc,install("step.fn"))); nprotect++;
    deltat = *(REAL(AS_NUMERIC(GET_SLOT(rproc,install("delta.t")))));
    PROTECT(X = euler_model_simulator(fn,xstart,times,params,deltat,method,
        zeronames,tcovar,covar,args,gnsi)); nprotect++;
  }
    break;
  case 5: // Gillespie's method
  {
    SEXP fn, vmatrix, hmax;
    PROTECT(fn = GET_SLOT(rproc,install("rate.fn"))); nprotect++;
    PROTECT(vmatrix = GET_SLOT(rproc,install("v"))); nprotect++;
    PROTECT(hmax = GET_SLOT(rproc,install("hmax"))); nprotect++;
    PROTECT(X = SSA_simulator(fn,xstart,times,params,vmatrix,tcovar,covar,
        zeronames,hmax,args,gnsi)); nprotect++;
  }
    break;
  case 1: default:
    errorcall(R_NilValue,"in 'rprocess': 'rprocess' is undefined.");
    break;
  }

  PROTECT(dimX = GET_DIM(X)); nprotect++;
  if ((isNull(dimX)) || (length(dimX) != 3)) {
    errorcall(R_NilValue,"in 'rprocess': user 'rprocess' must return a rank-3 array");
  }
  xdim = INTEGER(dimX);
  if ((xdim[0] != nvars) || (xdim[1] != nreps) || (xdim[2] != ntimes)) {
    errorcall(R_NilValue,"in 'rprocess': user 'rprocess' must return a %d x %d x %d array",nvars,nreps,ntimes);
  }
  if (isNull(GET_ROWNAMES(GET_DIMNAMES(X)))) {
    errorcall(R_NilValue,"in 'rprocess': user 'rprocess' must return an array with rownames");
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
