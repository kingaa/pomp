// -*- C++ -*-

#include <Rdefines.h>
#include <R_ext/Constants.h>
#include <string.h>

#include "pomp_internal.h"

static R_INLINE SEXP ret_array (int nvars, int nreps, int ntimes, SEXP Snames)
{
  SEXP X;
  int dim[3] = {nvars, nreps, ntimes};
  const char *dimnms[3] = {"variable","rep","time"};
  PROTECT(X = makearray(3,dim));
  setrownames(X,Snames,3);
  fixdimnames(X,dimnms,3);
  UNPROTECT(1);
  return X;
}

SEXP iterate_map (SEXP object, SEXP times, SEXP t0, SEXP x0, SEXP params, SEXP gnsi)
{
  int nprotect = 0;
  pompfunmode mode = undef;
  SEXP fn, args;
  SEXP X;
  SEXP Snames, Pnames, Cnames;
  SEXP skel, pompfun;
  SEXP zeronames;
  int *zidx = 0;
  int nvars, npars, nreps, nrepp, ntimes, ncovars, nzeros;
  int *dim;
  lookup_table_t covariate_table;
  double deltat, t;

  PROTECT(skel = GET_SLOT(object,install("skeleton"))); nprotect++;
  deltat = *(REAL(GET_SLOT(skel,install("delta.t"))));
  t = *(REAL(AS_NUMERIC(t0)));

  PROTECT(x0 = duplicate(x0)); nprotect++;
  PROTECT(x0 = as_matrix(x0)); nprotect++;
  dim = INTEGER(GET_DIM(x0));
  nvars = dim[0]; nreps = dim[1];

  PROTECT(params = as_matrix(params)); nprotect++;
  dim = INTEGER(GET_DIM(params));
  npars = dim[0]; nrepp = dim[1];

  PROTECT(times = AS_NUMERIC(times)); nprotect++;
  ntimes = LENGTH(times);

  PROTECT(Snames = GET_ROWNAMES(GET_DIMNAMES(x0))); nprotect++;
  PROTECT(Pnames = GET_ROWNAMES(GET_DIMNAMES(params))); nprotect++;
  PROTECT(Cnames = get_covariate_names(GET_SLOT(object,install("covar")))); nprotect++;

  // set up the covariate table
  covariate_table = make_covariate_table(GET_SLOT(object,install("covar")),&ncovars);

  // extract user-defined function
  PROTECT(pompfun = GET_SLOT(skel,install("skel.fn"))); nprotect++;
  PROTECT(fn = pomp_fun_handler(pompfun,gnsi,&mode,Snames,Pnames,NA_STRING,Cnames)); nprotect++;

  // extract 'userdata' as pairlist
  PROTECT(args = VectorToPairList(GET_SLOT(object,install("userdata")))); nprotect++;

  // get names and indices of accumulator variables
  PROTECT(zeronames = GET_SLOT(object,install("zeronames"))); nprotect++;
  nzeros = LENGTH(zeronames);
  if (nzeros > 0) {
    zidx = INTEGER(PROTECT(matchnames(Snames,zeronames,"state variables"))); nprotect++;
  }

  // create array to store results
  PROTECT(X = ret_array(nvars,nreps,ntimes,Snames)); nprotect++;

  // set up the computations
  switch (mode) {

  case Rfun: {

    PROTECT(args = add_skel_args(args,Snames,Pnames,Cnames)); nprotect++;

    iterate_skeleton_R(REAL(X),t,deltat,REAL(times),REAL(x0),REAL(params),
      fn,args,Snames,nvars,npars,ncovars,ntimes,nrepp,nreps,nzeros,
      &covariate_table,zidx);

  }

    break;

  case native: case regNative: {
    int *sidx, *pidx, *cidx;
    pomp_skeleton *ff;

    *((void **) (&ff)) = R_ExternalPtrAddr(fn);

    // construct state, parameter, covariate indices
    sidx = INTEGER(GET_SLOT(pompfun,install("stateindex")));
    pidx = INTEGER(GET_SLOT(pompfun,install("paramindex")));
    cidx = INTEGER(GET_SLOT(pompfun,install("covarindex")));

    iterate_skeleton_native(REAL(X),t,deltat,REAL(times),REAL(x0),REAL(params),
      nvars,npars,ncovars,ntimes,nrepp,nreps,nzeros,sidx,pidx,cidx,
      &covariate_table,zidx,ff,args);

  }

    break;

  default:

    break; // #nocov

  }

  UNPROTECT(nprotect);
  return X;
}

static struct {
  struct {
    pompfunmode mode;
    SEXP object;
    SEXP params;
    lookup_table_t covar_table;
    int nvars;
    int npars;
    int ncovars;
    int nreps;
  } common;
  union {
    struct {
      SEXP fn;
      SEXP args;
      SEXP Snames;
    } R_fun;
    struct {
      SEXP args;
      SEXP sindex;
      SEXP pindex;
      SEXP cindex;
      pomp_skeleton *fun;
    } native_code;
  } shared;
} _pomp_vf_eval_block;

#define COMMON(X) (_pomp_vf_eval_block.common.X)
#define RFUN(X)   (_pomp_vf_eval_block.shared.R_fun.X)
#define NAT(X)    (_pomp_vf_eval_block.shared.native_code.X)

SEXP pomp_desolve_setup (SEXP object, SEXP x0, SEXP params, SEXP gnsi) {
  int nprotect = 0;
  pompfunmode mode = undef;
  SEXP fn, args;
  SEXP Snames, Pnames, Cnames;
  SEXP pompfun;
  int *dim;
  int nvars, npars, nreps, ncovars;

  // extract user-defined skeleton function
  PROTECT(pompfun = GET_SLOT(GET_SLOT(object,install("skeleton")),install("skel.fn"))); nprotect++;
  // extract 'userdata' as pairlist
  PROTECT(args = VectorToPairList(GET_SLOT(object,install("userdata")))); nprotect++;

  COMMON(object) = object;
  COMMON(params) = params;
  if (!isNull(COMMON(object))) R_ReleaseObject(COMMON(object));
  if (!isNull(COMMON(params))) R_ReleaseObject(COMMON(params));
  R_PreserveObject(COMMON(object));
  R_PreserveObject(COMMON(params));

  dim = INTEGER(GET_DIM(x0));
  nvars = dim[0];
  COMMON(nvars) = nvars;

  dim = INTEGER(GET_DIM(params));
  npars = dim[0]; nreps = dim[1];
  COMMON(npars) = npars;
  COMMON(nreps) = nreps;

  // set up the covariate table
  COMMON(covar_table) = make_covariate_table(GET_SLOT(object,install("covar")),&ncovars);
  COMMON(ncovars) = ncovars;

  PROTECT(Snames = GET_ROWNAMES(GET_DIMNAMES(x0))); nprotect++;
  PROTECT(Pnames = GET_ROWNAMES(GET_DIMNAMES(params))); nprotect++;
  PROTECT(Cnames = get_covariate_names(GET_SLOT(object,install("covar")))); nprotect++;

  PROTECT(fn = pomp_fun_handler(pompfun,gnsi,&mode,Snames,Pnames,NA_STRING,Cnames)); nprotect++;

  COMMON(mode) = mode;

  switch (COMMON(mode)) {

  case Rfun: {

    PROTECT(RFUN(fn) = fn); nprotect++;
    PROTECT(RFUN(args) = add_skel_args(args,Snames,Pnames,Cnames)); nprotect++;
    PROTECT(RFUN(Snames) = Snames); nprotect++;

    if (!isNull(RFUN(fn))) R_ReleaseObject(RFUN(fn));
    if (!isNull(RFUN(args))) R_ReleaseObject(RFUN(args));
    if (!isNull(RFUN(Snames))) R_ReleaseObject(RFUN(Snames));

    R_PreserveObject(RFUN(fn));
    R_PreserveObject(RFUN(args));
    R_PreserveObject(RFUN(Snames));

  }

    break;

  case native: case regNative: {

    NAT(args) = args;

    PROTECT(NAT(sindex) = GET_SLOT(pompfun,install("stateindex"))); nprotect++;
    PROTECT(NAT(pindex) = GET_SLOT(pompfun,install("paramindex"))); nprotect++;
    PROTECT(NAT(cindex) = GET_SLOT(pompfun,install("covarindex"))); nprotect++;

    *((void **) (&(NAT(fun)))) = R_ExternalPtrAddr(fn);

    if (!isNull(NAT(args))) R_ReleaseObject(NAT(args));
    if (!isNull(NAT(sindex))) R_ReleaseObject(NAT(sindex));
    if (!isNull(NAT(pindex))) R_ReleaseObject(NAT(pindex));
    if (!isNull(NAT(cindex))) R_ReleaseObject(NAT(cindex));

    R_PreserveObject(NAT(args));
    R_PreserveObject(NAT(sindex));
    R_PreserveObject(NAT(pindex));
    R_PreserveObject(NAT(cindex));

  }

    break;

  default:

    errorcall(R_NilValue,"in 'pomp_desolve_setup': unrecognized 'mode'"); // # nocov

  }

  UNPROTECT(nprotect);

  return R_NilValue;
}

void pomp_vf_eval (int *neq, double *t, double *y, double *ydot, double *yout, int *ip)
{
  switch (COMMON(mode)) {

  case Rfun:			// R function

    eval_skeleton_R(ydot,t,y,REAL(COMMON(params)),
      RFUN(fn),RFUN(args),RFUN(Snames),
      COMMON(nvars),COMMON(npars),COMMON(ncovars),1,
      COMMON(nreps),COMMON(nreps),COMMON(nreps),
      &COMMON(covar_table));

    break;

  case native: case regNative:		// native code
    eval_skeleton_native(ydot,t,y,REAL(COMMON(params)),
      COMMON(nvars),COMMON(npars),COMMON(ncovars),1,
      COMMON(nreps),COMMON(nreps),COMMON(nreps),
      INTEGER(NAT(sindex)),INTEGER(NAT(pindex)),INTEGER(NAT(cindex)),
      &COMMON(covar_table),NAT(fun),NAT(args));

    break;

  default:

    errorcall(R_NilValue,"in 'pomp_vf_eval': unrecognized 'mode'"); // # nocov

  break;

  }
}

void pomp_desolve_takedown (void) {
  R_ReleaseObject(COMMON(object));
  R_ReleaseObject(COMMON(params));
  COMMON(object) = R_NilValue;
  COMMON(params) = R_NilValue;
  COMMON(nvars) = 0;
  COMMON(npars) = 0;
  COMMON(ncovars) = 0;
  COMMON(nreps) = 0;

  switch (COMMON(mode)) {

  case Rfun: {

    R_ReleaseObject(RFUN(fn));
    R_ReleaseObject(RFUN(args));
    R_ReleaseObject(RFUN(Snames));
    RFUN(fn) = R_NilValue;
    RFUN(args) = R_NilValue;
    RFUN(Snames) = R_NilValue;

  }

    break;

  case native: case regNative: {

    NAT(fun) = 0;
    R_ReleaseObject(NAT(args));
    R_ReleaseObject(NAT(sindex));
    R_ReleaseObject(NAT(pindex));
    R_ReleaseObject(NAT(cindex));
    NAT(args) = R_NilValue;
    NAT(sindex) = R_NilValue;
    NAT(pindex) = R_NilValue;
    NAT(cindex) = R_NilValue;

  }

    break;

  default:

    errorcall(R_NilValue,"in 'pomp_desolve_takedown': unrecognized 'mode'"); // # nocov

    break;

  }

  COMMON(mode) = undef;

}

#undef COMMON
#undef RFUN
#undef NAT
