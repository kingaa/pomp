// dear emacs, please treat this as -*- C++ -*-

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

#include "pomp_internal.h"

static SEXP pomp_default_rinit(SEXP params, SEXP Pnames,
  int npar, int nrep, int nsim);

SEXP do_rinit (SEXP object, SEXP params, SEXP t0, SEXP nsim, SEXP gnsi)
{
  int nprotect = 0;
  SEXP Pnames, Snames, pcnames;
  SEXP x = R_NilValue;
  int *dim;
  int npar, nrep, nvar, nsims, ns;
  int xdim[2];
  const char *dimnms[2] = {"variable","rep"};

  nsims = *(INTEGER(AS_INTEGER(nsim)));
  PROTECT(params = as_matrix(params)); nprotect++;
  PROTECT(Pnames = GET_ROWNAMES(GET_DIMNAMES(params))); nprotect++;
  PROTECT(pcnames = GET_COLNAMES(GET_DIMNAMES(params))); nprotect++;
  dim = INTEGER(GET_DIM(params));
  npar = dim[0]; nrep = dim[1];
  ns = nsims*nrep;

  SEXP pompfun, fcall, fn, covars = R_NilValue;
  pompfunmode mode = undef;
  double *cp = NULL;

  // extract the rinit function and its environment
  PROTECT(pompfun = GET_SLOT(object,install("rinit"))); nprotect++;
  PROTECT(fn = pomp_fun_handler(pompfun,gnsi,&mode)); nprotect++;

  // extract covariates and interpolate
  PROTECT(covars = lookup_in_table(GET_SLOT(object,install("covar")),t0)); nprotect++;
  if (LENGTH(covars) > 0) cp = REAL(covars);

  // extract userdata
  PROTECT(fcall = VectorToPairList(GET_SLOT(object,install("userdata")))); nprotect++;

  switch (mode) {
  case Rfun:			// use R function

  {
    SEXP par, rho, x1, x2;
    double *p, *pp, *xp, *xt;
    int j, *midx;

    // extract covariates and interpolate
    if (cp != NULL) { // add covars to call
      PROTECT(fcall = LCONS(covars,fcall)); nprotect++;
      SET_TAG(fcall,install("covars"));
    }

    // parameter vector
    PROTECT(par = NEW_NUMERIC(npar)); nprotect++;
    SET_NAMES(par,Pnames);
    pp = REAL(par);

    // finish constructing the call
    PROTECT(fcall = LCONS(t0,fcall)); nprotect++;
    SET_TAG(fcall,install("t0"));
    PROTECT(fcall = LCONS(par,fcall)); nprotect++;
    SET_TAG(fcall,install("params"));
    PROTECT(fcall = LCONS(fn,fcall)); nprotect++;

    // evaluation environment
    PROTECT(rho = (CLOENV(fn))); nprotect++;

    p = REAL(params);
    memcpy(pp,p,npar*sizeof(double));	   // copy the parameters
    PROTECT(x1 = eval(fcall,rho)); nprotect++; // do the call
    PROTECT(Snames = GET_NAMES(x1)); nprotect++;

    if (!IS_NUMERIC(x1) || isNull(Snames)) {
      errorcall(R_NilValue,"in 'rinit': user 'rinit' must return a named numeric vector");
    }

    nvar = LENGTH(x1);
    xp = REAL(x1);
    midx = INTEGER(PROTECT(match(Pnames,Snames,0))); nprotect++;

    for (j = 0; j < nvar; j++) {
      if (midx[j]!=0) {
        errorcall(R_NilValue,"in 'rinit': a state variable and a parameter share a single name: '%s'",CHAR(STRING_ELT(Snames,j)));
      }
    }

    xdim[0] = nvar; xdim[1] = ns;
    PROTECT(x = makearray(2,xdim)); nprotect++;
    setrownames(x,Snames,2);
    fixdimnames(x,dimnms,2);
    xt = REAL(x);

    memcpy(xt,xp,nvar*sizeof(double));

    for (j = 1, xt += nvar; j < ns; j++, xt += nvar) {
      memcpy(pp,p+npar*(j%nrep),npar*sizeof(double));
      PROTECT(x2 = eval(fcall,rho));
      xp = REAL(x2);
      if (LENGTH(x2)!=nvar)
        errorcall(R_NilValue,"in 'rinit': user 'rinit' returns vectors of non-uniform length");
      memcpy(xt,xp,nvar*sizeof(double));
      UNPROTECT(1);
    }

  }

    break;

  case native:		// use native routine

  {

    SEXP Cnames;
    int *sidx, *pidx, *cidx;
    double *xt, *ps, time;
    pomp_rinit *ff = NULL;
    int j;

    PROTECT(Snames = GET_SLOT(pompfun,install("statenames"))); nprotect++;
    PROTECT(Cnames = get_covariate_names(GET_SLOT(object,install("covar")))); nprotect++;

    // construct state, parameter, covariate, observable indices
    sidx = INTEGER(PROTECT(name_index(Snames,pompfun,"statenames","state variables"))); nprotect++;
    pidx = INTEGER(PROTECT(name_index(Pnames,pompfun,"paramnames","parameters"))); nprotect++;
    cidx = INTEGER(PROTECT(name_index(Cnames,pompfun,"covarnames","covariates"))); nprotect++;

    // address of native routine
    *((void **) (&ff)) = R_ExternalPtrAddr(fn);

    nvar = LENGTH(Snames);
    xdim[0] = nvar; xdim[1] = ns;
    PROTECT(x = makearray(2,xdim)); nprotect++;
    setrownames(x,Snames,2);
    fixdimnames(x,dimnms,2);

    set_pomp_userdata(fcall);
    GetRNGstate();

    time = *(REAL(t0));

    // loop over replicates
    for (j = 0, xt = REAL(x), ps = REAL(params); j < ns; j++, xt += nvar)
      (*ff)(xt,ps+npar*(j%nrep),time,sidx,pidx,cidx,cp);

    PutRNGstate();
    unset_pomp_userdata();

  }

    break;

  default:

    PROTECT(x = pomp_default_rinit(params,Pnames,npar,nrep,ns)); nprotect++;

    break;

  }

  if (nrep > 1) {
    SEXP c, dn, xn;
    int k, *p;

    if (isNull(pcnames)) {
      PROTECT(pcnames = NEW_INTEGER(nrep)); nprotect++;
      for (k = 0, p = INTEGER(pcnames); k < nrep; k++, p++) *p = k+1;
    }

    if (nsims > 1) {
      int k, *sp;

      PROTECT(xn = NEW_INTEGER(ns)); nprotect++;
      for (k = 0, sp = INTEGER(xn); k < ns; k++, sp++) *sp = (k/nrep)+1;
      PROTECT(c = LCONS(mkString("_"),R_NilValue)); nprotect++;
      SET_TAG(c,install("sep"));
      PROTECT(c = LCONS(install("paste"),LCONS(pcnames,LCONS(xn,c)))); nprotect++;
      PROTECT(xn = eval(c,R_BaseEnv)); nprotect++;
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
  int xdim[2];
  const char *dimnms[2] = {"variable","rep"};
  double *xp, *pp;
  int nprotect = 0;

  PROTECT(pat = NEW_CHARACTER(1)); nprotect++;
  SET_STRING_ELT(pat,0,mkChar("[\\_\\.]0$"));
  PROTECT(repl = NEW_CHARACTER(1)); nprotect++;
  SET_STRING_ELT(repl,0,mkChar(""));
  PROTECT(val = NEW_LOGICAL(1)); nprotect++;
  *(INTEGER(val)) = 1;

  // extract names of IVPs
  PROTECT(fcall = LCONS(val,R_NilValue)); nprotect++;
  SET_TAG(fcall,install("value"));
  PROTECT(fcall = LCONS(Pnames,fcall)); nprotect++;
  SET_TAG(fcall,install("x"));
  PROTECT(fcall = LCONS(pat,fcall)); nprotect++;
  SET_TAG(fcall,install("pattern"));
  PROTECT(fcall = LCONS(install("grep"),fcall)); nprotect++;
  PROTECT(ivpnames = eval(fcall,R_BaseEnv)); nprotect++;

  nvar = LENGTH(ivpnames);
  if (nvar < 1) {
    errorcall(R_NilValue,"in default 'rinit': there are no parameters with suffix '.0' or '_0'. See '?rinit_spec'.");
  }
  pidx = INTEGER(PROTECT(match(Pnames,ivpnames,0))); nprotect++;
  for (k = 0; k < nvar; k++) pidx[k]--;

  // construct names of state variables
  PROTECT(fcall = LCONS(ivpnames,R_NilValue)); nprotect++;
  SET_TAG(fcall,install("x"));
  PROTECT(fcall = LCONS(repl,fcall)); nprotect++;
  SET_TAG(fcall,install("replacement"));
  PROTECT(fcall = LCONS(pat,fcall)); nprotect++;
  SET_TAG(fcall,install("pattern"));
  PROTECT(fcall = LCONS(install("sub"),fcall)); nprotect++;
  PROTECT(statenames = eval(fcall,R_BaseEnv)); nprotect++;

  xdim[0] = nvar; xdim[1] = nsim;
  PROTECT(x = makearray(2,xdim)); nprotect++;
  setrownames(x,statenames,2);
  fixdimnames(x,dimnms,2);

  for (j = 0, xp = REAL(x); j < nsim; j++) {
    pp = REAL(params) + npar*(j%nrep);
    for (k = 0; k < nvar; k++, xp++) *xp = pp[pidx[k]];
  }

  UNPROTECT(nprotect);
  return x;
}
