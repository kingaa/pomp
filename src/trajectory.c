// -*- C++ -*-

#include <Rdefines.h>
#include <R_ext/Constants.h>

#include "pomp_internal.h"

void iterate_map_native (double *X, double *time, double *p,
			 double deltat, double t, double *x,
			 int ntimes, int nvars, int npars, int ncovars, int nzeros, int nreps,
			 int *sidx, int *pidx, int *cidx, int *zidx,
			 lookup_table *covar_table,
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
      table_lookup(covar_table,t,covars,0);
      for (j = 0, Xs = X, xs = x, ps = p; j < nreps; j++, Xs += nvars, xs += nvars, ps += npars) {
	(*ff)(Xs,xs,ps,sidx,pidx,cidx,ncovars,covars,t);
      }
      memcpy(x,X,nvars*nreps*sizeof(double));
      t += deltat;
    }
    if (nsteps == 0) memcpy(X,x,nvars*nreps*sizeof(double));
  }
  if (ncovars > 0) Free(covars);
  unset_pomp_userdata();
}

void iterate_map_R (double *X, double *time, double *p, 
		    double deltat, double t, double *x,
		    double *tp, double *xp, double *pp, double *cp,
		    int ntimes, int nvars, int npars, int nzeros, int nreps,
		    lookup_table *covar_table, int *zidx,
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
      table_lookup(covar_table,t,cp,0);
      for (j = 0, xs = x, ps = p; j < nreps; j++, xs += nvars, ps += npars) {
	*tp = t;
	memcpy(xp,xs,nvars*sizeof(double));
	memcpy(pp,ps,npars*sizeof(double));
	if (first) {
	  PROTECT(ans = eval(fcall,rho)); nprotect++;
	  if (LENGTH(ans)!=nvars)
	    errorcall(R_NilValue,"in 'trajectory': user 'skeleton' returns a vector of %d state variables but %d are expected",LENGTH(ans),nvars);
	  // get name information to fix possible alignment problems
	  PROTECT(nm = GET_NAMES(ans)); nprotect++;
	  use_names = !isNull(nm);
	  if (use_names) {
	    posn = INTEGER(PROTECT(matchnames(Snames,nm,"state variables"))); nprotect++;
	  }
	  fs = REAL(AS_NUMERIC(ans));
	  first = 0;
	} else {
	  fs = REAL(AS_NUMERIC(eval(fcall,rho)));
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
  SEXP pompfun;
  SEXP zeronames;
  int *zidx = 0;
  int nvars, npars, nreps, ntimes, ncovars, nzeros;
  int *dim;
  lookup_table covariate_table;
  double deltat, t;

  deltat = *(REAL(GET_SLOT(object,install("skelmap.delta.t"))));
  t = *(REAL(AS_NUMERIC(t0)));

  PROTECT(x0 = as_matrix(duplicate(x0))); nprotect++;
  dim = INTEGER(GET_DIM(x0));
  nvars = dim[0]; nreps = dim[1];

  PROTECT(params = as_matrix(params)); nprotect++;
  dim = INTEGER(GET_DIM(params));
  npars = dim[0];

  if (nreps != dim[1])
    errorcall(R_NilValue,"in 'trajectory': dimension mismatch between 'x0' and 'params'");

  PROTECT(times = AS_NUMERIC(times)); nprotect++;
  ntimes = LENGTH(times);

  PROTECT(Snames = GET_ROWNAMES(GET_DIMNAMES(x0))); nprotect++;
  PROTECT(Pnames = GET_ROWNAMES(GET_DIMNAMES(params))); nprotect++;
  PROTECT(Cnames = GET_COLNAMES(GET_DIMNAMES(GET_SLOT(object,install("covar"))))); nprotect++;
    
  // set up the covariate table
  covariate_table = make_covariate_table(object,&ncovars);

  // extract user-defined function
  PROTECT(pompfun = GET_SLOT(object,install("skeleton"))); nprotect++;
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
      int nprotect = 0;
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

      UNPROTECT(nprotect);
    }
    break;
  case native:                       // native skeleton
    {
      int nprotect = 0;
      int *sidx, *pidx, *cidx;
      pomp_skeleton *ff = (pomp_skeleton *) R_ExternalPtrAddr(fn);
      // construct state, parameter, covariate indices
      sidx = INTEGER(PROTECT(name_index(Snames,pompfun,"statenames"))); nprotect++;
      pidx = INTEGER(PROTECT(name_index(Pnames,pompfun,"paramnames"))); nprotect++;
      cidx = INTEGER(PROTECT(name_index(Cnames,pompfun,"covarnames"))); nprotect++;

      iterate_map_native(REAL(X),REAL(times),REAL(params),deltat,t,REAL(x0),
			 ntimes,nvars,npars,ncovars,nzeros,nreps,
			 sidx,pidx,cidx,zidx,&covariate_table,ff,args);

      UNPROTECT(nprotect);
    }
    break;
  default:
    errorcall(R_NilValue,"in 'iterate_map': unrecognized 'mode'");
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
    lookup_table covar_table;
    int nvars;
    int npars;
    int ncovars;
    int nreps;
  } common;
  union {
    struct {
      SEXP fcall;
      SEXP rho;
      SEXP Snames;
      SEXP tvec, xvec, pvec, cvec;
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
  PROTECT(pompfun = GET_SLOT(object,install("skeleton"))); nprotect++;
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
  COMMON(covar_table) = make_covariate_table(object,&ncovars);
  COMMON(ncovars) = ncovars;

  PROTECT(Snames = GET_ROWNAMES(GET_DIMNAMES(x0))); nprotect++;
  PROTECT(Pnames = GET_ROWNAMES(GET_DIMNAMES(params))); nprotect++;
  PROTECT(Cnames = GET_COLNAMES(GET_DIMNAMES(GET_SLOT(object,install("covar"))))); nprotect++;

  switch (COMMON(mode)) {
  case Rfun:			// R function
    // arguments of the R function
    PROTECT(RFUN(tvec) = NEW_NUMERIC(1)); nprotect++;
    PROTECT(RFUN(xvec) = NEW_NUMERIC(nvars)); nprotect++;
    PROTECT(RFUN(pvec) = NEW_NUMERIC(npars)); nprotect++;
    PROTECT(RFUN(cvec) = NEW_NUMERIC(ncovars)); nprotect++;
    SET_NAMES(RFUN(xvec),Snames);
    SET_NAMES(RFUN(pvec),Pnames);
    SET_NAMES(RFUN(cvec),Cnames);
    // set up the function call
    PROTECT(RFUN(fcall) = LCONS(RFUN(cvec),args)); nprotect++;
    SET_TAG(RFUN(fcall),install("covars"));
    PROTECT(RFUN(fcall) = LCONS(RFUN(pvec),RFUN(fcall))); nprotect++;
    SET_TAG(RFUN(fcall),install("params"));
    PROTECT(RFUN(fcall) = LCONS(RFUN(tvec),RFUN(fcall))); nprotect++;
    SET_TAG(RFUN(fcall),install("t"));
    PROTECT(RFUN(fcall) = LCONS(RFUN(xvec),RFUN(fcall))); nprotect++;
    SET_TAG(RFUN(fcall),install("x"));
    PROTECT(RFUN(fcall) = LCONS(fn,RFUN(fcall))); nprotect++;
    // environment of the user-defined function
    PROTECT(RFUN(rho) = (CLOENV(fn))); nprotect++;

    PROTECT(RFUN(Snames) = Snames); nprotect++;
    
    if (!isNull(RFUN(fcall))) R_ReleaseObject(RFUN(fcall));
    if (!isNull(RFUN(rho))) R_ReleaseObject(RFUN(rho));
    if (!isNull(RFUN(Snames))) R_ReleaseObject(RFUN(Snames));
    if (!isNull(RFUN(tvec))) R_ReleaseObject(RFUN(tvec));
    if (!isNull(RFUN(xvec))) R_ReleaseObject(RFUN(xvec));
    if (!isNull(RFUN(pvec))) R_ReleaseObject(RFUN(pvec));
    if (!isNull(RFUN(cvec))) R_ReleaseObject(RFUN(cvec));
    R_PreserveObject(RFUN(fcall));
    R_PreserveObject(RFUN(rho));
    R_PreserveObject(RFUN(Snames));
    R_PreserveObject(RFUN(tvec));
    R_PreserveObject(RFUN(xvec));
    R_PreserveObject(RFUN(pvec));
    R_PreserveObject(RFUN(cvec));

    break;
  case native:			// native code
    // set aside userdata
    NAT(args) = args;
    // construct index vectors
    PROTECT(NAT(sindex) = name_index(Snames,pompfun,"statenames")); nprotect++;
    PROTECT(NAT(pindex) = name_index(Pnames,pompfun,"paramnames")); nprotect++;
    PROTECT(NAT(cindex) = name_index(Cnames,pompfun,"covarnames")); nprotect++;
    // extract pointer to user-defined function
    NAT(fun) = (pomp_skeleton *) R_ExternalPtrAddr(fn);

    if (!isNull(NAT(args))) R_ReleaseObject(NAT(args));
    if (!isNull(NAT(sindex))) R_ReleaseObject(NAT(sindex));
    if (!isNull(NAT(pindex))) R_ReleaseObject(NAT(pindex));
    if (!isNull(NAT(cindex))) R_ReleaseObject(NAT(cindex));
    R_PreserveObject(NAT(args));
    R_PreserveObject(NAT(sindex));
    R_PreserveObject(NAT(pindex));
    R_PreserveObject(NAT(cindex));

    break;
  default:
    errorcall(R_NilValue,"in 'pomp_desolve_setup': unrecognized 'mode'");
    break;
  }
  UNPROTECT(nprotect);
  //  return retval;
  return R_NilValue;
}

void pomp_vf_eval (int *neq, double *t, double *y, double *ydot, double *yout, int *ip) 
{
  switch (COMMON(mode)) {
  case Rfun:			// R function
    eval_skeleton_R(ydot,t,y,REAL(COMMON(params)),
		    RFUN(fcall),RFUN(rho),RFUN(Snames),
		    REAL(RFUN(tvec)),REAL(RFUN(xvec)),REAL(RFUN(pvec)),REAL(RFUN(cvec)),
		    COMMON(nvars),COMMON(npars),1,COMMON(nreps),COMMON(nreps),COMMON(nreps),
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
    errorcall(R_NilValue,"in 'pomp_vf_eval': unrecognized 'mode'");
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
  case Rfun:			// R function
    R_ReleaseObject(RFUN(fcall));
    R_ReleaseObject(RFUN(rho));
    R_ReleaseObject(RFUN(Snames));
    R_ReleaseObject(RFUN(tvec));
    R_ReleaseObject(RFUN(xvec));
    R_ReleaseObject(RFUN(pvec));
    R_ReleaseObject(RFUN(cvec));
    RFUN(fcall) = R_NilValue;
    RFUN(rho) = R_NilValue;
    RFUN(Snames) = R_NilValue;
    RFUN(tvec) = R_NilValue;
    RFUN(xvec) = R_NilValue;
    RFUN(pvec) = R_NilValue;
    RFUN(cvec) = R_NilValue;
    break;
  case native:			// native code
    NAT(fun) = 0;
    R_ReleaseObject(NAT(args));
    R_ReleaseObject(NAT(sindex));
    R_ReleaseObject(NAT(pindex));
    R_ReleaseObject(NAT(cindex));
    NAT(args) = R_NilValue;
    NAT(sindex) = R_NilValue;
    NAT(pindex) = R_NilValue;
    NAT(cindex) = R_NilValue;
    break;
  default:
    errorcall(R_NilValue,"in 'pomp_desolve_takedown': unrecognized 'mode'");
    break;
  }
  COMMON(mode) = -1;
}

#undef COMMON
#undef RFUN
#undef NAT

// copy t(x[-1,-1]) -> y[,rep,]
// SEXP traj_transp_and_copy (SEXP y, SEXP x, SEXP rep) {
//   int nprotect = 0;
//   SEXP ans = R_NilValue;
//   int nvars, nreps, ntimes;
//   int i, j, k, m, n;
//   double *xp, *yp;

//   j = INTEGER(rep)[0]-1;
//   nvars = INTEGER(GET_DIM(y))[0];
//   nreps = INTEGER(GET_DIM(y))[1];
//   ntimes = INTEGER(GET_DIM(y))[2];
//   n = INTEGER(GET_DIM(x))[0];
//   m = nvars*nreps;

//   for (i = 0, xp = REAL(x)+n+1; i < nvars; i++, xp += n) {
//     for (k = 0, yp = REAL(y)+i+nvars*j; k < ntimes; k++, yp += m) {
//       *yp = xp[k];
//     }
//   }

//   UNPROTECT(nprotect);
//   return ans;
// }
