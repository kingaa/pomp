// dear emacs, please treat this as -*- C++ -*-

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

#include "pomp_internal.h"

// initializer.c
typedef double pomp_initializer(double *x, const double *p, double t, 
				const int *stateindex, const int *parindex, const int *covindex,
				const double *covars);

SEXP do_init_state (SEXP object, SEXP params, SEXP t0, SEXP gnsi)
{
  int nprotect = 0;
  SEXP Pnames, Snames, x;
  int *dim;
  int npar, nrep, nvar;
  int definit;
  int xdim[2];
  const char *dimnms[2] = {"variable","rep"};

  PROTECT(params = as_matrix(params)); nprotect++;
  PROTECT(Pnames = GET_ROWNAMES(GET_DIMNAMES(params))); nprotect++;
  dim = INTEGER(GET_DIM(params));
  npar = dim[0]; nrep = dim[1]; 

  definit = *(INTEGER(GET_SLOT(object,install("default.init"))));

  if (definit) {		// default initializer

    SEXP fcall, pat, repl, val, ivpnames, statenames;
    int *pidx, j, k;
    double *xp, *pp;
  
    PROTECT(pat = NEW_CHARACTER(1)); nprotect++;
    SET_STRING_ELT(pat,0,mkChar("\\.0$"));
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
      error("default initializer error: no parameter names ending in '.0' found: see 'pomp' documentation");
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

    xdim[0] = nvar; xdim[1] = nrep;
    PROTECT(x = makearray(2,xdim)); nprotect++;
    setrownames(x,statenames,2);
    fixdimnames(x,dimnms,2);

    for (j = 0, xp = REAL(x), pp = REAL(params); j < nrep; j++, pp += npar)
      for (k = 0; k < nvar; k++, xp++)
	*xp = pp[pidx[k]];

  } else {			// user-supplied initializer
    
    SEXP pompfun, fcall, fn, tcovar, covar, covars;
    pompfunmode mode = undef;
    double *cp = NULL;

    // extract the initializer function and its environment
    PROTECT(pompfun = GET_SLOT(object,install("initializer"))); nprotect++;
    PROTECT(fn = pomp_fun_handler(pompfun,gnsi,&mode)); nprotect++;
    
    // extract covariates and interpolate
    PROTECT(tcovar = GET_SLOT(object,install("tcovar"))); nprotect++;
    if (LENGTH(tcovar) > 0) {	// do table lookup
      PROTECT(covar = GET_SLOT(object,install("covar"))); nprotect++;
      PROTECT(covars = lookup_in_table(tcovar,covar,t0)); nprotect++;
      cp = REAL(covars);
    }
	
    // extract userdata
    PROTECT(fcall = VectorToPairList(GET_SLOT(object,install("userdata")))); nprotect++;
	
    switch (mode) {
    case Rfun:			// use R function

      {
	SEXP par, rho, x1, x2, mindex;
	double *p, *pp, *xp, *xt;
	int j, *midx;

	// extract covariates and interpolate
	if (LENGTH(tcovar) > 0) { // add covars to call
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
	  UNPROTECT(nprotect);
	  error("'init.state' error: user 'initializer' must return a named numeric vector");
	}
	
	nvar = LENGTH(x1);
	xp = REAL(x1);
	midx = INTEGER(PROTECT(match(Pnames,Snames,0))); nprotect++;
	
	for (j = 0; j < nvar; j++) {
	  if (midx[j]!=0) {
	    UNPROTECT(nprotect);
	    error("a state variable and a parameter share a single name: '%s'",CHARACTER_DATA(STRING_ELT(Snames,j)));
	  }
	}
	
	xdim[0] = nvar; xdim[1] = nrep;
	PROTECT(x = makearray(2,xdim)); nprotect++;
	setrownames(x,Snames,2);
	fixdimnames(x,dimnms,2);
	xt = REAL(x);
	
	memcpy(xt,xp,nvar*sizeof(double));
	
	for (j = 1, p += npar, xt += nvar; j < nrep; j++, p += npar, xt += nvar) {
	  memcpy(pp,p,npar*sizeof(double));
	  PROTECT(x2 = eval(fcall,rho));
	  xp = REAL(x2);
	  if (LENGTH(x2)!=nvar)
	    error("user initializer returns vectors of non-uniform length");
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
	pomp_initializer *ff = NULL;
	int j;

	PROTECT(Snames = GET_SLOT(pompfun,install("statenames"))); nprotect++;
	PROTECT(Cnames = GET_SLOT(pompfun,install("covarnames"))); nprotect++;

	// construct state, parameter, covariate, observable indices
	sidx = INTEGER(PROTECT(name_index(Snames,pompfun,"statenames"))); nprotect++;
	pidx = INTEGER(PROTECT(name_index(Pnames,pompfun,"paramnames"))); nprotect++;
	cidx = INTEGER(PROTECT(name_index(Cnames,pompfun,"covarnames"))); nprotect++;
	
	// address of native routine
	ff = (pomp_initializer *) R_ExternalPtrAddr(fn);
	
	nvar = LENGTH(Snames);
	xdim[0] = nvar; xdim[1] = nrep;
	PROTECT(x = makearray(2,xdim)); nprotect++;
	setrownames(x,Snames,2);
	fixdimnames(x,dimnms,2);
	
	set_pomp_userdata(fcall);
	GetRNGstate();

	time = *(REAL(t0));

	// loop over replicates
	for (j = 0, xt = REAL(x), ps = REAL(params); j < nrep; j++, xt += nvar, ps += npar)
	  (*ff)(xt,ps,time,sidx,pidx,cidx,cp);

	PutRNGstate();
	unset_pomp_userdata();
      
      }

      break;
      
    default:
      
      error("unrecognized 'mode' slot in 'init.state'");
      break;

    }

  }

  UNPROTECT(nprotect);
  return x;
}
