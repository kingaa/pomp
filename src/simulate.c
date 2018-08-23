// -*- C++ -*-

#include <Rdefines.h>
#include <string.h>

#include "pomp_internal.h"

SEXP do_simulate (SEXP object, SEXP params, SEXP nsim, SEXP gnsi)
{
  int nprotect = 0;
  SEXP t0, times, alltimes, xstart, x, y;
  SEXP ans, ans_names;
  SEXP statenames, paramnames, simnames, offset;
  int ntimes;

  PROTECT(offset = NEW_INTEGER(1)); nprotect++;
  *(INTEGER(offset)) = 1;

  if (LENGTH(nsim)<1)
    errorcall(R_NilValue,"'nsim' must be a single integer");
  if (LENGTH(nsim)>1)
    warningcall(R_NilValue,"in 'simulate': only the first number in 'nsim' is significant");

  PROTECT(params = as_matrix(params)); nprotect++;
  PROTECT(paramnames = GET_ROWNAMES(GET_DIMNAMES(params))); nprotect++;

  PROTECT(t0 = GET_SLOT(object,install("t0"))); nprotect++;
  PROTECT(times = GET_SLOT(object,install("times"))); nprotect++;
  ntimes = LENGTH(times);
  PROTECT(alltimes = NEW_NUMERIC(ntimes+1)); nprotect++;
  *(REAL(alltimes)) = *(REAL(t0));			// copy t0 into alltimes[1]
  memcpy(REAL(alltimes)+1,REAL(times),ntimes*sizeof(double)); // copy times into alltimes[-1]

  // initialize the simulations
  PROTECT(xstart = do_rinit(object,params,t0,nsim,gnsi)); nprotect++;
  PROTECT(statenames = GET_ROWNAMES(GET_DIMNAMES(xstart))); nprotect++;
  PROTECT(simnames = GET_COLNAMES(GET_DIMNAMES(xstart))); nprotect++;

  // call 'rprocess' to simulate state process
  PROTECT(x = do_rprocess(object,xstart,alltimes,params,offset,gnsi)); nprotect++;
  // call 'rmeasure' to simulate the measurement process
  PROTECT(y = do_rmeasure(object,x,times,params,gnsi)); nprotect++;
  setcolnames(x,simnames);
  setcolnames(y,simnames);

  PROTECT(ans = NEW_LIST(2)); nprotect++;
  PROTECT(ans_names = NEW_CHARACTER(2)); nprotect++;
  SET_STRING_ELT(ans_names,0,mkChar("states"));
  SET_STRING_ELT(ans_names,1,mkChar("obs"));
  SET_NAMES(ans,ans_names);
  SET_ELEMENT(ans,0,x);
  SET_ELEMENT(ans,1,y);
  UNPROTECT(nprotect);
  return ans;

}
