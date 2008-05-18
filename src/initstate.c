// dear emacs, please treat this as -*- C++ -*-

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

SEXP do_init_state (SEXP object, SEXP params, SEXP t0)
{
  int nprotect = 0;
  SEXP X, fn, slotname, userdata, fcall, rho;
  PROTECT(slotname = NEW_CHARACTER(1)); nprotect++;
  // extract the process function
  SET_STRING_ELT(slotname,0,mkChar("initializer"));
  PROTECT(fn = GET_SLOT(object,slotname)); nprotect++;
  // extract the userdata
  SET_STRING_ELT(slotname,0,mkChar("userdata"));
  PROTECT(userdata = GET_SLOT(object,slotname)); nprotect++;
  // construct the call
  PROTECT(fcall = LCONS(fn,LCONS(params,LCONS(t0,VectorToPairList(userdata))))); nprotect++;
  PROTECT(rho = (CLOENV(fn))); nprotect++; // environment of the function
  PROTECT(X = eval(fcall,rho)); nprotect++; // do the call
  if (!IS_NUMERIC(X) || isNull(GET_NAMES(X)))
    error("init.state error: user 'initializer' must return a named numeric vector");
  UNPROTECT(nprotect);
  return X;
}
