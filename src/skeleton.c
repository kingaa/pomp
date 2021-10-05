// dear emacs, please treat this as -*- C++ -*-

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Arith.h>

#include "pomp_internal.h"

SEXP add_skel_args(SEXP args, SEXP Snames, SEXP Pnames, SEXP Cnames);
void eval_skeleton_R(double *f, double *time, double *x, double *p, SEXP fn, SEXP args, SEXP Snames, int nvars, int npars, int ncovars, int ntimes, int nrepx, int nrepp, int nreps, lookup_table_t *covar_table, double *cov);
void iterate_skeleton_R(double *X, double t, double deltat, double *time, double *x, double *p, SEXP fn, SEXP args, SEXP Snames, int nvars, int npars, int ncovars, int ntimes, int nrepp, int nreps, int nzeros, lookup_table_t *covar_table, int *zeroindex, double *cov);
void eval_skeleton_native(double *f, double *time, double *x, double *p, int nvars, int npars, int ncovars, int ntimes, int nrepx, int nrepp, int nreps, int *sidx, int *pidx, int *cidx, lookup_table_t *covar_table, pomp_skeleton *fun, SEXP args, double *cov);
void iterate_skeleton_native(double *X, double t, double deltat, double *time, double *x, double *p, int nvars, int npars, int ncovars, int ntimes, int nrepp, int nreps, int nzeros, int *sidx, int *pidx, int *cidx, lookup_table_t *covar_table, int *zeroindex, pomp_skeleton *fun, SEXP args, double *cov);
SEXP do_skeleton(SEXP object, SEXP x, SEXP t, SEXP params, SEXP gnsi);

SEXP add_skel_args (SEXP args, SEXP Snames, SEXP Pnames, SEXP Cnames)
{

  SEXP var;
  int v;

  PROTECT(args);

  // Covariates
  for (v = LENGTH(Cnames)-1; v >= 0; v--) {
    var = NEW_NUMERIC(1);
    args = LCONS(var,args);
    UNPROTECT(1);
    PROTECT(args);
    SET_TAG(args,installChar(STRING_ELT(Cnames,v)));
  }

  // Parameters
  for (v = LENGTH(Pnames)-1; v >= 0; v--) {
    var = NEW_NUMERIC(1);
    args = LCONS(var,args);
    UNPROTECT(1);
    PROTECT(args);
    SET_TAG(args,installChar(STRING_ELT(Pnames,v)));
  }

  // Latent state variables
  for (v = LENGTH(Snames)-1; v >= 0; v--) {
    var = NEW_NUMERIC(1);
    args = LCONS(var,args);
    UNPROTECT(1);
    PROTECT(args);
    SET_TAG(args,installChar(STRING_ELT(Snames,v)));
  }

  // Time
  var = NEW_NUMERIC(1);
  args = LCONS(var,args);
  UNPROTECT(1);
  PROTECT(args);
  SET_TAG(args,install("t"));

  UNPROTECT(1);
  return args;

}

static R_INLINE SEXP eval_call (
    SEXP fn, SEXP args,
    double *t,
    double *x, int nvar,
    double *p, int npar,
    double *c, int ncov)
{

  SEXP var = args, ans, ob;
  int v;

  *(REAL(CAR(var))) = *t; var = CDR(var);
  for (v = 0; v < nvar; v++, x++, var=CDR(var)) *(REAL(CAR(var))) = *x;
  for (v = 0; v < npar; v++, p++, var=CDR(var)) *(REAL(CAR(var))) = *p;
  for (v = 0; v < ncov; v++, c++, var=CDR(var)) *(REAL(CAR(var))) = *c;

  PROTECT(ob = LCONS(fn,args));
  PROTECT(ans = eval(ob,CLOENV(fn)));

  UNPROTECT(2);
  return ans;

}

static R_INLINE SEXP ret_array (int nvars, int nreps, int ntimes, SEXP Snames)
{
  SEXP X;
  int dim[3] = {nvars, nreps, ntimes};
  const char *dimnms[3] = {"variable",".id","time"};
  PROTECT(X = makearray(3,dim));
  setrownames(X,Snames,3);
  fixdimnames(X,dimnms,3);
  UNPROTECT(1);
  return X;
}

void eval_skeleton_R (
    double *f, double *time, double *x, double *p,
    SEXP fn, SEXP args, SEXP Snames,
    int nvars, int npars, int ncovars, int ntimes,
    int nrepx, int nrepp, int nreps,
    lookup_table_t *covar_table, double *cov)
{

  SEXP ans, nm;
  double *fs;
  int *posn = 0;
  int i, j, k;
  int first = 1;
  int nprotect = 0;
  
  for (k = 0; k < ntimes; k++, time++) {

    R_CheckUserInterrupt();

    // interpolate the covar functions for the covariates
    table_lookup(covar_table,*time,cov);

    for (j = 0; j < nreps; j++, f += nvars) {

      if (first) {

        PROTECT(ans = eval_call(fn,args,time,
          x+nvars*((j%nrepx)+nrepx*k),nvars,
          p+npars*(j%nrepp),npars,
          cov,ncovars));

          if (LENGTH(ans)!=nvars)
            err("'skeleton' returns a vector of %d state variables but %d are expected.",LENGTH(ans),nvars);

          // get name information to fix alignment problems
          PROTECT(nm = GET_NAMES(ans));
          if (invalid_names(nm))
            err("'skeleton' must return a named numeric vector.");
          posn = INTEGER(PROTECT(matchnames(Snames,nm,"state variables")));

          fs = REAL(AS_NUMERIC(ans));
          for (i = 0; i < nvars; i++) f[posn[i]] = fs[i];

	  nprotect += 3;
	  first = 0;

      } else {

        PROTECT(ans = eval_call(fn,args,time,
          x+nvars*((j%nrepx)+nrepx*k),nvars,
          p+npars*(j%nrepp),npars,
          cov,ncovars));

        fs = REAL(AS_NUMERIC(ans));
        for (i = 0; i < nvars; i++) f[posn[i]] = fs[i];

        UNPROTECT(1);

      }

    }
  }

  UNPROTECT(nprotect);

}

void iterate_skeleton_R (
    double *X, double t, double deltat,
    double *time, double *x, double *p,
    SEXP fn, SEXP args, SEXP Snames,
    int nvars, int npars, int ncovars, int ntimes,
    int nrepp, int nreps, int nzeros,
    lookup_table_t *covar_table, int *zeroindex,
    double *cov)
{

  SEXP ans, nm;
  double *ap, *xs;
  int nsteps;
  int *posn = 0;
  int h, i, j, k;
  int first = 1;
  int nprotect = 0;

  for (k = 0; k < ntimes; k++, time++, X += nvars*nreps) {

    R_CheckUserInterrupt();

    // compute number of steps needed
    nsteps = num_map_steps(t,*time,deltat);

    // set accumulator variables to zero
    for (i = 0; i < nzeros; i++)
      for (j = 0, xs = &x[zeroindex[i]]; j < nreps; j++, xs += nvars)
        *xs = 0.0;

    for (h = 0; h < nsteps; h++) {

      // interpolate the covariates
      table_lookup(covar_table,t,cov);

      for (j = 0, xs = x; j < nreps; j++, xs += nvars) {

        if (first) {

          PROTECT(ans = eval_call(fn,args,&t,xs,nvars,p+npars*(j%nrepp),npars,cov,ncovars));

          if (LENGTH(ans) != nvars)
            err("'skeleton' returns a vector of %d state variables but %d are expected.",LENGTH(ans),nvars);

          // get name information to fix alignment problems
          PROTECT(nm = GET_NAMES(ans));
          if (invalid_names(nm)) err("'skeleton' must return a named numeric vector.");
          posn = INTEGER(PROTECT(matchnames(Snames,nm,"state variables")));

          ap = REAL(AS_NUMERIC(ans));
          for (i = 0; i < nvars; i++) xs[posn[i]] = ap[i];

	  nprotect += 3;
          first = 0;

        } else {

          PROTECT(ans = eval_call(fn,args,&t,xs,nvars,p+npars*(j%nrepp),npars,cov,ncovars));

          ap = REAL(AS_NUMERIC(ans));
          for (i = 0; i < nvars; i++) xs[posn[i]] = ap[i];

          UNPROTECT(1);

        }

      }

      if (h != nsteps-1) {
        t += deltat;
      } else {
        deltat = *time - t;
        t = *time;
      }

    }

    memcpy(X,x,nvars*nreps*sizeof(double));

  }

  UNPROTECT(nprotect);

}

void eval_skeleton_native (
    double *f, double *time, double *x, double *p,
    int nvars, int npars, int ncovars, int ntimes,
    int nrepx, int nrepp, int nreps,
    int *sidx, int *pidx, int *cidx,
    lookup_table_t *covar_table, pomp_skeleton *fun, SEXP args,
    double *cov)
{
  double *xp, *pp;
  int j, k;

  set_pomp_userdata(args);

  for (k = 0; k < ntimes; k++, time++) {

    R_CheckUserInterrupt();

    table_lookup(covar_table,*time,cov);

    for (j = 0; j < nreps; j++, f += nvars) {

      xp = &x[nvars*((j%nrepx)+nrepx*k)];
      pp = &p[npars*(j%nrepp)];

      (*fun)(f,xp,pp,sidx,pidx,cidx,cov,*time);

    }
  }

  unset_pomp_userdata();

}

void iterate_skeleton_native (
    double *X, double t, double deltat,
    double *time, double *x, double *p,
    int nvars, int npars, int ncovars, int ntimes,
    int nrepp, int nreps, int nzeros,
    int *sidx, int *pidx, int *cidx,
    lookup_table_t *covar_table, int *zeroindex,
    pomp_skeleton *fun, SEXP args, double *cov)
{
  int nsteps;
  double *xs, *Xs;
  int h, i, j, k;

  set_pomp_userdata(args);

  for (k = 0; k < ntimes; k++, time++, X += nvars*nreps) {

    R_CheckUserInterrupt();

    // compute number of steps needed
    nsteps = num_map_steps(t,*time,deltat);

    // set accumulator variables to zero
    for (i = 0; i < nzeros; i++)
      for (j = 0, xs = &x[zeroindex[i]]; j < nreps; j++, xs += nvars)
        *xs = 0.0;

    for (h = 0; h < nsteps; h++) {

      // interpolate the covariates
      table_lookup(covar_table,t,cov);

      for (j = 0, Xs = X, xs = x; j < nreps; j++, Xs += nvars, xs += nvars) {

        (*fun)(Xs,xs,p+npars*(j%nrepp),sidx,pidx,cidx,cov,t);

      }

      memcpy(x,X,nvars*nreps*sizeof(double));

      if (h != nsteps-1) {
        t += deltat;
      } else {
        deltat = *time - t;
        t = *time;
      }

    }

    if (nsteps == 0) memcpy(X,x,nvars*nreps*sizeof(double));

  }

  unset_pomp_userdata();

}

SEXP do_skeleton (SEXP object, SEXP x, SEXP t, SEXP params, SEXP gnsi)
{

  int nvars, npars, nrepp, nrepx, nreps, ntimes, ncovars;
  pompfunmode mode = undef;
  int *dim;
  SEXP Snames, Cnames, Pnames;
  SEXP pompfun, cov, ob;
  SEXP fn, args, F;
  lookup_table_t covariate_table;

  PROTECT(t = AS_NUMERIC(t));
  ntimes = LENGTH(t);

  PROTECT(x = as_state_array(x));
  dim = INTEGER(GET_DIM(x));
  nvars = dim[0]; nrepx = dim[1];
  if (ntimes != dim[2])
    err("length of 't' and 3rd dimension of 'x' do not agree.");

  PROTECT(params = as_matrix(params));
  dim = INTEGER(GET_DIM(params));
  npars = dim[0]; nrepp = dim[1];

  // 2nd dimension of 'x' and 'params' need not agree
  nreps = (nrepp > nrepx) ? nrepp : nrepx;
  if ((nreps % nrepp != 0) || (nreps % nrepx != 0))
    err("2nd dimensions of 'x' and 'params' are incompatible");

  PROTECT(Snames = GET_ROWNAMES(GET_DIMNAMES(x)));
  PROTECT(Pnames = GET_ROWNAMES(GET_DIMNAMES(params)));
  PROTECT(Cnames = get_covariate_names(GET_SLOT(object,install("covar"))));

  // set up the covariate table
  covariate_table = make_covariate_table(GET_SLOT(object,install("covar")),&ncovars);
  PROTECT(cov = NEW_NUMERIC(ncovars));

  // extract the user-defined function
  PROTECT(ob = GET_SLOT(object,install("skeleton")));
  PROTECT(pompfun = GET_SLOT(ob,install("skel.fn")));
  PROTECT(fn = pomp_fun_handler(pompfun,gnsi,&mode,Snames,Pnames,NA_STRING,Cnames));

  // extract 'userdata' as pairlist
  PROTECT(args = VectorToPairList(GET_SLOT(object,install("userdata"))));

  PROTECT(F = ret_array(nvars,nreps,ntimes,Snames));

  int nprotect = 12;
  
  switch (mode) {

  case Rfun: {

    PROTECT(args = add_skel_args(args,Snames,Pnames,Cnames)); nprotect++;

    eval_skeleton_R(REAL(F),REAL(t),REAL(x),REAL(params),
      fn,args,Snames,
      nvars,npars,ncovars,ntimes,nrepx,nrepp,nreps,
      &covariate_table,REAL(cov));

  }

    break;

  case native: case regNative: {
    int *sidx, *pidx, *cidx;
    pomp_skeleton *ff = NULL;

    sidx = INTEGER(GET_SLOT(pompfun,install("stateindex")));
    pidx = INTEGER(GET_SLOT(pompfun,install("paramindex")));
    cidx = INTEGER(GET_SLOT(pompfun,install("covarindex")));

    *((void **) (&ff)) = R_ExternalPtrAddr(fn);

    eval_skeleton_native(
      REAL(F),REAL(t),REAL(x),REAL(params),
      nvars,npars,ncovars,ntimes,nrepx,nrepp,nreps,
      sidx,pidx,cidx,&covariate_table,ff,args,REAL(cov));

  }

    break;

  default: {

    double *ft = REAL(F);
    int i, n = nvars*nreps*ntimes;
    for (i = 0; i < n; i++, ft++) *ft = R_NaReal;
    warn("'skeleton' unspecified: NAs generated.");

  }

  }

  UNPROTECT(nprotect);
  return F;
}
