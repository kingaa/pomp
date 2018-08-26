// dear emacs, please treat this as -*- C++ -*-

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

#include "pomp_internal.h"

static SEXP pomp_default_rinit(SEXP params, SEXP Pnames,
  int npar, int nrep, int nsim);

static R_INLINE SEXP add_args (SEXP args, SEXP Pnames, SEXP Cnames)
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

  // Time
  PROTECT(var = NEW_NUMERIC(1)); nprotect++;
  PROTECT(args = LCONS(var,args)); nprotect++;
  SET_TAG(args,install("t0"));

  UNPROTECT(nprotect);
  return args;

}

static R_INLINE SEXP eval_call (SEXP fn, SEXP args,
  double *t0, double *p, int npar, double *c, int ncov)
{

  SEXP var = args, ans;
  int v;

  *(REAL(CAR(var))) = *t0; var = CDR(var);
  for (v = 0; v < npar; v++, p++, var=CDR(var)) *(REAL(CAR(var))) = *p;
  for (v = 0; v < ncov; v++, c++, var=CDR(var)) *(REAL(CAR(var))) = *c;

  PROTECT(ans = eval(LCONS(fn,args),CLOENV(fn)));

  UNPROTECT(1);
  return ans;

}

static R_INLINE SEXP ret_array (int m, int n, SEXP names)
{
  int dim[2] = {m, n};
  const char *dimnm[2] = {"variable","rep"};
  SEXP X;
  PROTECT(X = makearray(2,dim));
  setrownames(X,names,2);
  fixdimnames(X,dimnm,2);
  UNPROTECT(1);
  return X;
}

SEXP do_rinit (SEXP object, SEXP params, SEXP t0, SEXP nsim, SEXP gnsi)
{
  int nprotect = 0;
  SEXP Pnames, Cnames, Snames, pcnames;
  SEXP x = R_NilValue;
  SEXP pompfun, fn, args;
  pompfunmode mode = undef;
  lookup_table_t covariate_table;
  int *dim;
  int npar, nrep, nvar, ncovars, nsims, ns;

  nsims = *(INTEGER(AS_INTEGER(nsim)));
  PROTECT(params = as_matrix(params)); nprotect++;
  PROTECT(Pnames = GET_ROWNAMES(GET_DIMNAMES(params))); nprotect++;
  PROTECT(pcnames = GET_COLNAMES(GET_DIMNAMES(params))); nprotect++;

  dim = INTEGER(GET_DIM(params));
  npar = dim[0]; nrep = dim[1];
  ns = nsims*nrep;

  // extract the rinit function and its environment
  PROTECT(pompfun = GET_SLOT(object,install("rinit"))); nprotect++;
  PROTECT(fn = pomp_fun_handler(pompfun,gnsi,&mode)); nprotect++;

  // set up the covariate table
  covariate_table = make_covariate_table(GET_SLOT(object,install("covar")),&ncovars);
  PROTECT(Cnames = get_covariate_names(GET_SLOT(object,install("covar")))); nprotect++;

  double cov[ncovars];
  table_lookup(&covariate_table,*(REAL(t0)),cov);

  // extract userdata
  PROTECT(args = VectorToPairList(GET_SLOT(object,install("userdata")))); nprotect++;

  switch (mode) {
  case Rfun: {

    SEXP ans;
    double *time = REAL(AS_NUMERIC(t0));
    double *ps = REAL(params);
    double *xs, *xt = NULL;
    int *midx;
    int j;

    PROTECT(args = add_args(args,Pnames,Cnames)); nprotect++;

    PROTECT(ans = eval_call(fn,args,time,ps,npar,cov,ncovars)); nprotect++;
    PROTECT(Snames = GET_NAMES(ans)); nprotect++;

    if (!IS_NUMERIC(ans) || isNull(Snames))
      errorcall(R_NilValue,"user 'rinit' must return a named numeric vector");

    nvar = LENGTH(ans);
    xs = REAL(ans);
    midx = INTEGER(PROTECT(match(Pnames,Snames,0))); nprotect++;

    for (j = 0; j < nvar; j++) {
      if (midx[j] != 0)
        errorcall(R_NilValue,"a state variable and a parameter share the name: '%s'",CHAR(STRING_ELT(Snames,j)));
    }

    PROTECT(x = ret_array(nvar,ns,Snames)); nprotect++;
    xt = REAL(x);

    memcpy(xt,xs,nvar*sizeof(double));

    for (j = 1, xt += nvar, ps += npar; j < ns; j++, xt += nvar, ps += npar) {
      PROTECT(ans = eval_call(fn,args,time,ps,npar,cov,ncovars));
      xs = REAL(ans);
      if (LENGTH(ans) != nvar)
        errorcall(R_NilValue,"user 'rinit' returns vectors of non-uniform length");
      memcpy(xt,xs,nvar*sizeof(double));
      UNPROTECT(1);
    }

  }

    break;

  case native: case regNative: {

    int *sidx, *pidx, *cidx;
    double *xt, *ps, time;
    pomp_rinit *ff = NULL;
    int j;

    PROTECT(Snames = GET_SLOT(pompfun,install("statenames"))); nprotect++;
    nvar = LENGTH(Snames);
    PROTECT(x = ret_array(nvar,ns,Snames)); nprotect++;

    sidx = INTEGER(PROTECT(name_index(Snames,pompfun,"statenames","state variables"))); nprotect++;
    pidx = INTEGER(PROTECT(name_index(Pnames,pompfun,"paramnames","parameters"))); nprotect++;
    cidx = INTEGER(PROTECT(name_index(Cnames,pompfun,"covarnames","covariates"))); nprotect++;

    // address of native routine
    *((void **) (&ff)) = R_ExternalPtrAddr(fn);

    set_pomp_userdata(args);
    GetRNGstate();

    time = *(REAL(t0));

    // loop over replicates
    for (j = 0, xt = REAL(x), ps = REAL(params); j < ns; j++, xt += nvar)
      (*ff)(xt,ps+npar*(j%nrep),time,sidx,pidx,cidx,cov);

    PutRNGstate();
    unset_pomp_userdata();

  }

    break;

  default: {

    PROTECT(x = pomp_default_rinit(params,Pnames,npar,nrep,ns)); nprotect++;

  }

  break;

  }

  // now add column names
  if (nrep > 1) {
    SEXP dn, xn;
    int k, *p;

    if (isNull(pcnames)) {
      PROTECT(pcnames = NEW_INTEGER(nrep)); nprotect++;
      for (k = 0, p = INTEGER(pcnames); k < nrep; k++, p++) *p = k+1;
    }

    if (nsims > 1) {
      int k, *sp;

      PROTECT(xn = NEW_INTEGER(ns)); nprotect++;
      for (k = 0, sp = INTEGER(xn); k < ns; k++, sp++) *sp = (k/nrep)+1;
      PROTECT(xn = paste(pcnames,xn,mkString("_"))); nprotect++;
      PROTECT(dn = GET_DIMNAMES(x)); nprotect++;
      SET_ELEMENT(dn,1,xn);
      SET_DIMNAMES(x,dn);

    } else {

      PROTECT(dn = GET_DIMNAMES(x)); nprotect++;
      SET_ELEMENT(dn,1,pcnames);
      SET_DIMNAMES(x,dn);

    }
  }

  UNPROTECT(nprotect);
  return x;
}

static SEXP pomp_default_rinit (SEXP params, SEXP Pnames,
  int npar, int nrep, int nsim)
{
  SEXP fcall, pat, repl, val, ivpnames, statenames, x;
  int *pidx;
  int nvar, j, k;
  double *xp, *pp;
  int nprotect = 0;

  // set up search pattern: ".0" or "_0"
  PROTECT(pat = NEW_CHARACTER(1)); nprotect++;
  SET_STRING_ELT(pat,0,mkChar("[\\_\\.]0$"));
  PROTECT(repl = NEW_CHARACTER(1)); nprotect++;
  SET_STRING_ELT(repl,0,mkChar(""));
  PROTECT(val = NEW_LOGICAL(1)); nprotect++;
  *(INTEGER(val)) = 1;

  // extract names of IVPs using 'grep'
  PROTECT(fcall = LCONS(val,R_NilValue)); nprotect++;
  SET_TAG(fcall,install("value"));
  PROTECT(fcall = LCONS(Pnames,fcall)); nprotect++;
  SET_TAG(fcall,install("x"));
  PROTECT(fcall = LCONS(pat,fcall)); nprotect++;
  SET_TAG(fcall,install("pattern"));
  PROTECT(fcall = LCONS(install("grep"),fcall)); nprotect++;
  PROTECT(ivpnames = eval(fcall,R_BaseEnv)); nprotect++;

  nvar = LENGTH(ivpnames);
  if (nvar < 1)
    errorcall(R_NilValue,"in default 'rinit': there are no parameters with suffix '.0' or '_0'. See '?rinit_spec'.");

  pidx = INTEGER(PROTECT(match(Pnames,ivpnames,0))); nprotect++;
  for (k = 0; k < nvar; k++) pidx[k]--;

  // construct names of state variables using 'sub'
  PROTECT(fcall = LCONS(ivpnames,R_NilValue)); nprotect++;
  SET_TAG(fcall,install("x"));
  PROTECT(fcall = LCONS(repl,fcall)); nprotect++;
  SET_TAG(fcall,install("replacement"));
  PROTECT(fcall = LCONS(pat,fcall)); nprotect++;
  SET_TAG(fcall,install("pattern"));
  PROTECT(fcall = LCONS(install("sub"),fcall)); nprotect++;
  PROTECT(statenames = eval(fcall,R_BaseEnv)); nprotect++;

  PROTECT(x = ret_array(nvar,nsim,statenames)); nprotect++;

  for (j = 0, xp = REAL(x); j < nsim; j++) {
    pp = REAL(params) + npar*(j%nrep);
    for (k = 0; k < nvar; k++, xp++) *xp = pp[pidx[k]];
  }

  UNPROTECT(nprotect);
  return x;
}
