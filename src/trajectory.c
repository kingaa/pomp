// -*- C++ -*-

#include <Rdefines.h>
#include "internal.h"

static R_INLINE SEXP ret_array (int nvars, int nreps, int ntimes, SEXP Snames, SEXP repnames)
{
  SEXP X;
  int dim[3] = {nvars, nreps, ntimes};
  const char *dimnms[3] = {"name",".id","time"};
  PROTECT(X = makearray(3,dim));
  setrownames(X,Snames,3);
  setcolnames(X,repnames);
  fixdimnames(X,dimnms,3);
  UNPROTECT(1);
  return X;
}

SEXP iterate_map (SEXP object, SEXP times, SEXP t0, SEXP x0, SEXP params, SEXP gnsi)
{

  pompfunmode mode = undef;
  SEXP fn, args;
  SEXP X, cov;
  SEXP Snames, Pnames, Cnames;
  SEXP skel, pompfun;
  SEXP accumvars;
  SEXP repnames;
  int *zidx = 0;
  int nvars, npars, nreps, nrepp, ntimes, ncovars, nzeros;
  int *dim;
  lookup_table_t covariate_table;
  double deltat, t;

  PROTECT(skel = GET_SLOT(object,install("skeleton")));
  deltat = *(REAL(GET_SLOT(skel,install("delta.t"))));
  t = *(REAL(AS_NUMERIC(t0)));

  PROTECT(x0 = duplicate(x0));
  PROTECT(x0 = as_matrix(x0));
  dim = INTEGER(GET_DIM(x0));
  nvars = dim[0]; nreps = dim[1];

  PROTECT(params = as_matrix(params));
  dim = INTEGER(GET_DIM(params));
  npars = dim[0]; nrepp = dim[1];
  PROTECT(repnames = GET_COLNAMES(GET_DIMNAMES(params)));

  PROTECT(times = AS_NUMERIC(times));
  ntimes = LENGTH(times);

  PROTECT(Snames = GET_ROWNAMES(GET_DIMNAMES(x0)));
  PROTECT(Pnames = GET_ROWNAMES(GET_DIMNAMES(params)));
  PROTECT(Cnames = get_covariate_names(GET_SLOT(object,install("covar"))));

  // set up the covariate table
  covariate_table = make_covariate_table(GET_SLOT(object,install("covar")),&ncovars);
  PROTECT(cov = NEW_NUMERIC(ncovars));

  // extract user-defined function
  PROTECT(pompfun = GET_SLOT(skel,install("skel.fn")));
  PROTECT(fn = pomp_fun_handler(pompfun,gnsi,&mode,Snames,Pnames,NA_STRING,Cnames));

  // extract 'userdata' as pairlist
  PROTECT(args = GET_SLOT(object,install("userdata")));

  // create array to store results
  PROTECT(X = ret_array(nvars,nreps,ntimes,Snames,repnames));

  // get names and indices of accumulator variables
  PROTECT(accumvars = GET_SLOT(object,install("accumvars")));
  nzeros = LENGTH(accumvars);

  int nprotect = 15;

  if (nzeros > 0) {
    zidx = INTEGER(PROTECT(matchnames(Snames,accumvars,"state variables"))); nprotect++;
  }

  // set up the computations
  switch (mode) {

  case Rfun: {

    PROTECT(args = add_skel_args(args,Snames,Pnames,Cnames)); nprotect++;

    iterate_skeleton_R(REAL(X),t,deltat,REAL(times),REAL(x0),REAL(params),
                       fn,args,Snames,nvars,npars,ncovars,ntimes,nrepp,nreps,nzeros,
                       &covariate_table,zidx,REAL(cov));

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
                            &covariate_table,zidx,ff,args,REAL(cov));

  }

    break;

  default: // #nocov

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
    SEXP cov;
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

  pompfunmode mode = undef;
  SEXP fn, args;
  SEXP Snames, Pnames, Cnames;
  SEXP pompfun, ob;
  int *dim;
  int nvars, npars, nreps, ncovars;

  // extract user-defined skeleton function
  PROTECT(ob = GET_SLOT(object,install("skeleton")));
  PROTECT(pompfun = GET_SLOT(ob,install("skel.fn")));
  // extract 'userdata' as pairlist
  PROTECT(args = GET_SLOT(object,install("userdata")));

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
  PROTECT(COMMON(cov) = NEW_NUMERIC(ncovars));
  R_PreserveObject(COMMON(cov));

  PROTECT(Snames = GET_ROWNAMES(GET_DIMNAMES(x0)));
  PROTECT(Pnames = GET_ROWNAMES(GET_DIMNAMES(params)));
  PROTECT(Cnames = get_covariate_names(GET_SLOT(object,install("covar"))));

  PROTECT(fn = pomp_fun_handler(pompfun,gnsi,&mode,Snames,Pnames,NA_STRING,Cnames));

  COMMON(mode) = mode;

  int nprotect = 8;

  switch (COMMON(mode)) {

  case Rfun: {

    PROTECT(RFUN(fn) = fn);
    PROTECT(RFUN(args) = add_skel_args(args,Snames,Pnames,Cnames));
    PROTECT(RFUN(Snames) = Snames);
    nprotect += 3;

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

    PROTECT(NAT(sindex) = GET_SLOT(pompfun,install("stateindex")));
    PROTECT(NAT(pindex) = GET_SLOT(pompfun,install("paramindex")));
    PROTECT(NAT(cindex) = GET_SLOT(pompfun,install("covarindex")));
    nprotect += 3;

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

  default: // #nocov

    err("in 'pomp_desolve_setup': unrecognized 'mode'"); // # nocov

  }

  UNPROTECT(nprotect);

  return R_NilValue;
}

void pomp_vf_eval (int *neq, double *t, double *y, double *ydot, double *yout, int *ip)
{
  switch (COMMON(mode)) {

  case Rfun:                    // R function

    eval_skeleton_R(ydot,t,y,REAL(COMMON(params)),
                    RFUN(fn),RFUN(args),RFUN(Snames),
                    COMMON(nvars),COMMON(npars),COMMON(ncovars),1,
                    COMMON(nreps),COMMON(nreps),COMMON(nreps),
                    &COMMON(covar_table),REAL(COMMON(cov)));

    break;

  case native: case regNative:          // native code
    eval_skeleton_native(ydot,t,y,REAL(COMMON(params)),
                         COMMON(nvars),COMMON(npars),COMMON(ncovars),1,
                         COMMON(nreps),COMMON(nreps),COMMON(nreps),
                         INTEGER(NAT(sindex)),INTEGER(NAT(pindex)),INTEGER(NAT(cindex)),
                         &COMMON(covar_table),NAT(fun),NAT(args),REAL(COMMON(cov)));

    break;

  default: // #nocov

    err("in 'pomp_vf_eval': unrecognized 'mode'"); // # nocov

    break;

  }
}

SEXP pomp_desolve_takedown (void) {
  R_ReleaseObject(COMMON(object));
  R_ReleaseObject(COMMON(params));
  R_ReleaseObject(COMMON(cov));
  COMMON(object) = R_NilValue;
  COMMON(params) = R_NilValue;
  COMMON(cov) = R_NilValue;
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

  default: // #nocov

    err("in 'pomp_desolve_takedown': unrecognized 'mode'"); // # nocov

    break;

  }

  COMMON(mode) = undef;

  return R_NilValue;
}

#undef COMMON
#undef RFUN
#undef NAT
