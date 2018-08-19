// dear emacs, please treat this as -*- C++ -*-

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "pomp_internal.h"

void _flat_improper_dprior (double *lik, double *p, int give_log, int *parindex) {
  *lik = (give_log) ? 0.0 : 1.0;
}

SEXP do_dprior (SEXP object, SEXP params, SEXP log, SEXP gnsi)
{
  int nprotect = 0;
  pompfunmode mode = undef;
  int npars, nreps;
  SEXP Pnames, F, fn, fcall;
  SEXP pompfun;
  int *dim;

  PROTECT(params = as_matrix(params)); nprotect++;
  dim = INTEGER(GET_DIM(params));
  npars = dim[0]; nreps = dim[1];

  PROTECT(Pnames = GET_ROWNAMES(GET_DIMNAMES(params))); nprotect++;

  // extract the user-defined function
  PROTECT(pompfun = GET_SLOT(object,install("dprior"))); nprotect++;
  PROTECT(fn = pomp_fun_handler(pompfun,gnsi,&mode)); nprotect++;

  // extract 'userdata' as pairlist
  PROTECT(fcall = VectorToPairList(GET_SLOT(object,install("userdata")))); nprotect++;

  // to store results
  PROTECT(F = NEW_NUMERIC(nreps)); nprotect++;

  // first do setup
  switch (mode) {
  case Rfun:			// use R function

  {
    SEXP pvec, rho;
    double *pp, *ps, *pt;
    int j;

    // temporary storage
    PROTECT(pvec = NEW_NUMERIC(npars)); nprotect++;
    SET_NAMES(pvec,Pnames);

    // set up the function call
    PROTECT(fcall = LCONS(AS_LOGICAL(log),fcall)); nprotect++;
    SET_TAG(fcall,install("log"));
    PROTECT(fcall = LCONS(pvec,fcall)); nprotect++;
    SET_TAG(fcall,install("params"));
    PROTECT(fcall = LCONS(fn,fcall)); nprotect++;

    // get the function's environment
    PROTECT(rho = (CLOENV(fn))); nprotect++;

    pp = REAL(pvec);

    for (j = 0, ps = REAL(params), pt = REAL(F); j < nreps; j++, ps += npars, pt++) {

      memcpy(pp,ps,npars*sizeof(double));

      *pt = *(REAL(AS_NUMERIC(PROTECT(eval(fcall,rho)))));
      UNPROTECT(1);

    }
  }

    break;

  case native:			// use native routine

  {
    int give_log, *pidx = 0;
    pomp_dprior *ff = NULL;
    double *ps, *pt;
    int j;

    // construct state, parameter, covariate, observable indices
    pidx = INTEGER(PROTECT(name_index(Pnames,pompfun,"paramnames","parameters"))); nprotect++;

    // address of native routine
    *((void **) (&ff)) = R_ExternalPtrAddr(fn);

    give_log = *(INTEGER(AS_INTEGER(log)));

    R_CheckUserInterrupt();	// check for user interrupt

    set_pomp_userdata(fcall);

    // loop over replicates
    for (j = 0, pt = REAL(F), ps = REAL(params); j < nreps; j++, ps += npars, pt++)
      (*ff)(pt,ps,give_log,pidx);

    unset_pomp_userdata();
  }

    break;

  default:

    errorcall(R_NilValue,"in 'dprior': unrecognized 'mode'"); // # nocov

  break;

  }

  UNPROTECT(nprotect);
  return F;
}
