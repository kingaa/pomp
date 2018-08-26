// -*- C++ -*-

#include <Rdefines.h>
#include <R_ext/Constants.h>

#include "pomp_internal.h"

static void iterate_map_native (double *X, double *time, double *p,
  double deltat, double t, double *x,
  int ntimes, int nvars, int npars, int ncovars, int nzeros, int nreps,
  int *sidx, int *pidx, int *cidx, int *zidx,
  lookup_table_t *covar_table,
  pomp_skeleton *ff, SEXP args)
{
  double *covars = NULL;
  int nsteps;
  double *Xs, *xs, *ps;
  int h, i, j, k;
  set_pomp_userdata(args);
  if (ncovars > 0) covars = (double *) Calloc(ncovars, double);
  for (k = 0; k < ntimes; k++, time++, X += nvars*nreps) {
    R_CheckUserInterrupt();
    for (i = 0; i < nzeros; i++)
      for (j = 0, xs = &x[zidx[i]]; j < nreps; j++, xs += nvars)
        *xs = 0.0;
    nsteps = num_map_steps(t,*time,deltat);
    for (h = 0; h < nsteps; h++) {
      table_lookup(covar_table,t,covars);
      for (j = 0, Xs = X, xs = x, ps = p; j < nreps; j++, Xs += nvars, xs += nvars, ps += npars) {
        (*ff)(Xs,xs,ps,sidx,pidx,cidx,covars,t);
      }
      memcpy(x,X,nvars*nreps*sizeof(double));
      t += deltat;
    }
    if (nsteps == 0) memcpy(X,x,nvars*nreps*sizeof(double));
  }
  if (ncovars > 0) Free(covars);
  unset_pomp_userdata();
}

static void iterate_map_R (double *X, double *time, double *p,
  double deltat, double t, double *x,
  double *tp, double *xp, double *pp, double *cp,
  int ntimes, int nvars, int npars, int nzeros, int nreps,
  lookup_table_t *covar_table, int *zidx,
  SEXP Snames, SEXP fcall, SEXP rho)
{
  int nprotect = 0;
  int first = 1;
  int use_names;
  int nsteps;
  SEXP ans, nm;
  double *fs, *xs, *ps;
  int *posn = 0;
  int h, i, j, k;

  for (k = 0; k < ntimes; k++, time++, X += nvars*nreps) {
    R_CheckUserInterrupt();
    nsteps = num_map_steps(t,*time,deltat);
    for (i = 0; i < nzeros; i++)
      for (j = 0, xs = &x[zidx[i]]; j < nreps; j++, xs += nvars)
        *xs = 0.0;
    for (h = 0; h < nsteps; h++) {
      table_lookup(covar_table,t,cp);
      for (j = 0, xs = x, ps = p; j < nreps; j++, xs += nvars, ps += npars) {
        *tp = t;
        memcpy(xp,xs,nvars*sizeof(double));
        memcpy(pp,ps,npars*sizeof(double));
        if (first) {
          PROTECT(ans = eval(fcall,rho)); nprotect++;
          if (LENGTH(ans)!=nvars)
            errorcall(R_NilValue,"user 'skeleton' returns a vector of %d state variables but %d are expected",LENGTH(ans),nvars);
          // get name information to fix possible alignment problems
          PROTECT(nm = GET_NAMES(ans)); nprotect++;
          use_names = !isNull(nm);
          if (use_names) {
            posn = INTEGER(PROTECT(matchnames(Snames,nm,"state variables"))); nprotect++;
          }
          fs = REAL(AS_NUMERIC(ans));
          first = 0;
        } else {
          fs = REAL(AS_NUMERIC(PROTECT(eval(fcall,rho))));
          UNPROTECT(1);
        }
        if (use_names)
          for (i = 0; i < nvars; i++) xs[posn[i]] = fs[i];
        else
          for (i = 0; i < nvars; i++) xs[i] = fs[i];
      }
      t += deltat;
    }
    memcpy(X,x,nvars*nreps*sizeof(double));
  }
  UNPROTECT(nprotect);
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
  int nvars, npars, nreps, ntimes, ncovars, nzeros;
  int *dim;
  lookup_table_t covariate_table;
  double deltat, t;

  PROTECT(skel = GET_SLOT(object,install("skeleton"))); nprotect++;
  deltat = *(REAL(GET_SLOT(skel,install("delta.t"))));
  t = *(REAL(AS_NUMERIC(t0)));

  PROTECT(x0 = as_matrix(PROTECT(duplicate(x0)))); nprotect += 2;
  dim = INTEGER(GET_DIM(x0));
  nvars = dim[0]; nreps = dim[1];

  PROTECT(params = as_matrix(params)); nprotect++;
  dim = INTEGER(GET_DIM(params));
  npars = dim[0];

  if (nreps != dim[1])
    errorcall(R_NilValue,"dimension mismatch between 'x0' and 'params'"); // # nocov

  PROTECT(times = AS_NUMERIC(times)); nprotect++;
  ntimes = LENGTH(times);

  PROTECT(Snames = GET_ROWNAMES(GET_DIMNAMES(x0))); nprotect++;
  PROTECT(Pnames = GET_ROWNAMES(GET_DIMNAMES(params))); nprotect++;
  PROTECT(Cnames = get_covariate_names(GET_SLOT(object,install("covar")))); nprotect++;

  // set up the covariate table
  covariate_table = make_covariate_table(GET_SLOT(object,install("covar")),&ncovars);

  // extract user-defined function
  PROTECT(pompfun = GET_SLOT(skel,install("skel.fn"))); nprotect++;
  PROTECT(fn = pomp_fun_handler(pompfun,gnsi,&mode)); nprotect++;

  // extract 'userdata' as pairlist
  PROTECT(args = VectorToPairList(GET_SLOT(object,install("userdata")))); nprotect++;

  // get names and indices of accumulator variables
  PROTECT(zeronames = GET_SLOT(object,install("zeronames"))); nprotect++;
  nzeros = LENGTH(zeronames);
  if (nzeros > 0) {
    zidx = INTEGER(PROTECT(matchnames(Snames,zeronames,"state variables"))); nprotect++;
  }

  // create array to store results
  {
    int dim[3] = {nvars, nreps, ntimes};
    PROTECT(X = makearray(3,dim)); nprotect++;
    setrownames(X,Snames,3);
  }

  // set up the computations
  switch (mode) {

  case Rfun:                       // R function
  {
    SEXP cvec, tvec, xvec, pvec;
    SEXP fcall, rho;
    PROTECT(tvec = NEW_NUMERIC(1)); nprotect++;
    PROTECT(xvec = NEW_NUMERIC(nvars)); nprotect++;
    PROTECT(pvec = NEW_NUMERIC(npars)); nprotect++;
    PROTECT(cvec = NEW_NUMERIC(ncovars)); nprotect++;
    SET_NAMES(xvec,Snames);
    SET_NAMES(pvec,Pnames);
    SET_NAMES(cvec,Cnames);
    // set up the function call
    PROTECT(fcall = LCONS(cvec,args)); nprotect++;
    SET_TAG(fcall,install("covars"));
    PROTECT(fcall = LCONS(pvec,fcall)); nprotect++;
    SET_TAG(fcall,install("params"));
    PROTECT(fcall = LCONS(tvec,fcall)); nprotect++;
    SET_TAG(fcall,install("t"));
    PROTECT(fcall = LCONS(xvec,fcall)); nprotect++;
    SET_TAG(fcall,install("x"));
    PROTECT(fcall = LCONS(fn,fcall)); nprotect++;
    // get function's environment
    PROTECT(rho = (CLOENV(fn))); nprotect++;

    iterate_map_R(REAL(X),REAL(times),REAL(params),
      deltat,t,REAL(x0),
      REAL(tvec),REAL(xvec),REAL(pvec),REAL(cvec),
      ntimes,nvars,npars,nzeros,nreps,
      &covariate_table,zidx,Snames,fcall,rho);
  }

    break;

  case native:                       // native skeleton
  {
    int *sidx, *pidx, *cidx;
    pomp_skeleton *ff;
    *((void **) (&ff)) = R_ExternalPtrAddr(fn);
    // construct state, parameter, covariate indices
    sidx = INTEGER(PROTECT(name_index(Snames,pompfun,"statenames","state variables"))); nprotect++;
    pidx = INTEGER(PROTECT(name_index(Pnames,pompfun,"paramnames","parameters"))); nprotect++;
    cidx = INTEGER(PROTECT(name_index(Cnames,pompfun,"covarnames","covariates"))); nprotect++;

    iterate_map_native(REAL(X),REAL(times),REAL(params),deltat,t,REAL(x0),
      ntimes,nvars,npars,ncovars,nzeros,nreps,
      sidx,pidx,cidx,zidx,&covariate_table,ff,args);

  }

    break;

  default:

    errorcall(R_NilValue,"unrecognized 'mode'"); // # nocov

  break;

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
  PROTECT(fn = pomp_fun_handler(pompfun,gnsi,&mode)); nprotect++;
  // extract 'userdata' as pairlist
  PROTECT(args = VectorToPairList(GET_SLOT(object,install("userdata")))); nprotect++;

  COMMON(mode) = mode;

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

    PROTECT(NAT(sindex) = name_index(Snames,pompfun,"statenames","state variables")); nprotect++;
    PROTECT(NAT(pindex) = name_index(Pnames,pompfun,"paramnames","parameters")); nprotect++;
    PROTECT(NAT(cindex) = name_index(Cnames,pompfun,"covarnames","covariates")); nprotect++;

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

  break;

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

  case native:			// native code
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
