// dear emacs, please treat this as -*- C++ -*-

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Arith.h>

#include "pomp_internal.h"

SEXP add_skel_args (SEXP args, SEXP Snames, SEXP Pnames, SEXP Cnames)
{
  int nprotect = 0;
  SEXP var;
  int v;

  // Covariates
  for (v = LENGTH(Cnames)-1; v >= 0; v--) {
    PROTECT(var = NEW_NUMERIC(1)); nprotect++;
    PROTECT(args = LCONS(var,args)); nprotect++;
    SET_TAG(args,install(CHAR(STRING_ELT(Cnames,v))));
  }

  // Parameters
  for (v = LENGTH(Pnames)-1; v >= 0; v--) {
    PROTECT(var = NEW_NUMERIC(1)); nprotect++;
    PROTECT(args = LCONS(var,args)); nprotect++;
    SET_TAG(args,install(CHAR(STRING_ELT(Pnames,v))));
  }

  // Latent state variables
  for (v = LENGTH(Snames)-1; v >= 0; v--) {
    PROTECT(var = NEW_NUMERIC(1)); nprotect++;
    PROTECT(args = LCONS(var,args)); nprotect++;
    SET_TAG(args,install(CHAR(STRING_ELT(Snames,v))));
  }

  // Time
  PROTECT(var = NEW_NUMERIC(1)); nprotect++;
  PROTECT(args = LCONS(var,args)); nprotect++;
  SET_TAG(args,install("t"));

  UNPROTECT(nprotect);
  return args;

}

static R_INLINE SEXP eval_call (
    SEXP fn, SEXP args,
    double *t,
    double *x, int nvar,
    double *p, int npar,
    double *c, int ncov)
{

  SEXP var = args, ans;
  int v;

  *(REAL(CAR(var))) = *t; var = CDR(var);
  for (v = 0; v < nvar; v++, x++, var=CDR(var)) *(REAL(CAR(var))) = *x;
  for (v = 0; v < npar; v++, p++, var=CDR(var)) *(REAL(CAR(var))) = *p;
  for (v = 0; v < ncov; v++, c++, var=CDR(var)) *(REAL(CAR(var))) = *c;

  PROTECT(ans = eval(LCONS(fn,args),CLOENV(fn)));

  UNPROTECT(1);
  return ans;

}

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

void eval_skeleton_R (
    double *f, double *time, double *x, double *p,
    SEXP fn, SEXP args, SEXP Snames,
    int nvars, int npars, int ncovars, int ntimes,
    int nrepx, int nrepp, int nreps,
    lookup_table_t *covar_table)
{
  int nprotect = 0;
  SEXP ans, nm;
  double *fs;
  int *posn = 0;
  double cov[ncovars];
  int i, j, k;

  for (k = 0; k < ntimes; k++, time++) {

    R_CheckUserInterrupt();

    // interpolate the covar functions for the covariates
    table_lookup(covar_table,*time,cov);

    for (j = 0; j < nreps; j++, f += nvars) {

      if (k == 0 && j == 0) {

        PROTECT(ans = eval_call(fn,args,time,
          x+nvars*((j%nrepx)+nrepx*k),nvars,
          p+npars*(j%nrepp),npars,
          cov,ncovars)); nprotect++;

        if (LENGTH(ans)!=nvars)
          errorcall(R_NilValue,"'skeleton' returns a vector of %d state variables but %d are expected.",LENGTH(ans),nvars);

        // get name information to fix alignment problems
        PROTECT(nm = GET_NAMES(ans)); nprotect++;
        if (isNull(nm))
          errorcall(R_NilValue,"'skeleton' must return a named numeric vector.");
        posn = INTEGER(PROTECT(matchnames(Snames,nm,"state variables"))); nprotect++;
        fs = REAL(AS_NUMERIC(ans));

        for (i = 0; i < nvars; i++) f[posn[i]] = fs[i];

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

void eval_skeleton_native (
    double *f, double *time, double *x, double *p,
    int nvars, int npars, int ncovars, int ntimes,
    int nrepx, int nrepp, int nreps,
    int *sidx, int *pidx, int *cidx,
    lookup_table_t *covar_table, pomp_skeleton *fun, SEXP args)
{
  double *xp, *pp;
  double cov[ncovars];
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

SEXP do_skeleton (SEXP object, SEXP x, SEXP t, SEXP params, SEXP gnsi)
{
  int nprotect = 0;
  int nvars, npars, nrepp, nrepx, nreps, ntimes, ncovars;
  pompfunmode mode = undef;
  int *dim;
  SEXP Snames, Cnames, Pnames;
  SEXP pompfun;
  SEXP fn, args, F;
  lookup_table_t covariate_table;

  PROTECT(t = AS_NUMERIC(t)); nprotect++;
  ntimes = LENGTH(t);

  PROTECT(x = as_state_array(x)); nprotect++;
  dim = INTEGER(GET_DIM(x));
  nvars = dim[0]; nrepx = dim[1];
  if (ntimes != dim[2])
    errorcall(R_NilValue,"length of 't' and 3rd dimension of 'x' do not agree");

  PROTECT(params = as_matrix(params)); nprotect++;
  dim = INTEGER(GET_DIM(params));
  npars = dim[0]; nrepp = dim[1];

  // 2nd dimension of 'x' and 'params' need not agree
  nreps = (nrepp > nrepx) ? nrepp : nrepx;
  if ((nreps % nrepp != 0) || (nreps % nrepx != 0))
    errorcall(R_NilValue,"2nd dimensions of 'x' and 'params' are incompatible");

  PROTECT(Snames = GET_ROWNAMES(GET_DIMNAMES(x))); nprotect++;
  PROTECT(Pnames = GET_ROWNAMES(GET_DIMNAMES(params))); nprotect++;
  PROTECT(Cnames = get_covariate_names(GET_SLOT(object,install("covar")))); nprotect++;

  // set up the covariate table
  covariate_table = make_covariate_table(GET_SLOT(object,install("covar")),&ncovars);

  // extract the user-defined function
  PROTECT(pompfun = GET_SLOT(GET_SLOT(object,install("skeleton")),install("skel.fn"))); nprotect++;
  PROTECT(fn = pomp_fun_handler(pompfun,gnsi,&mode)); nprotect++;

  // extract 'userdata' as pairlist
  PROTECT(args = VectorToPairList(GET_SLOT(object,install("userdata")))); nprotect++;

  PROTECT(F = ret_array(nvars,nreps,ntimes,Snames)); nprotect++;

  switch (mode) {

  case Rfun: {

    PROTECT(args = add_skel_args(args,Snames,Pnames,Cnames)); nprotect++;

    eval_skeleton_R(REAL(F),REAL(t),REAL(x),REAL(params),
      fn,args,Snames,
      nvars,npars,ncovars,ntimes,nrepx,nrepp,nreps,
      &covariate_table);

  }

    break;

  case native: {
    int *sidx, *pidx, *cidx;
    pomp_skeleton *ff = NULL;

    sidx = INTEGER(PROTECT(name_index(Snames,pompfun,"statenames","state variables"))); nprotect++;
    pidx = INTEGER(PROTECT(name_index(Pnames,pompfun,"paramnames","parameters"))); nprotect++;
    cidx = INTEGER(PROTECT(name_index(Cnames,pompfun,"covarnames","covariates"))); nprotect++;

    *((void **) (&ff)) = R_ExternalPtrAddr(fn);

    eval_skeleton_native(
      REAL(F),REAL(t),REAL(x),REAL(params),
      nvars,npars,ncovars,ntimes,nrepx,nrepp,nreps,
      sidx,pidx,cidx,&covariate_table,ff,args);

  }

    break;

  default: {

    double *ft = REAL(F);
    int i, n = nvars*nreps*ntimes;
    for (i = 0; i < n; i++, ft++) *ft = R_NaReal;
  }

  break;

  }

  UNPROTECT(nprotect);
  return F;
}
