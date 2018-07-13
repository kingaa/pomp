// dear emacs, please treat this as -*- C++ -*-

#include "pomp_internal.h"

static SEXP __pomp_ptr_userdata;
#define USERDATA  (__pomp_ptr_userdata)

void set_pomp_userdata (SEXP userdata) {
  USERDATA = userdata;
}

const SEXP get_pomp_userdata (const char *name) {
  SEXP elt = getPairListElement(USERDATA,name);
  if (isNull(elt)) errorcall(R_NilValue,"no user-provided element '%s' is found",name);
  return elt;
}

const int *get_pomp_userdata_int (const char *name) {
  SEXP elt = getPairListElement(USERDATA,name);
  if (isNull(elt)) errorcall(R_NilValue,"no user-provided element named '%s' is found",name);
  if (!isInteger(elt)) errorcall(R_NilValue,"user-provided element '%s' is not an integer",name);
  return INTEGER(elt);
}

const double *get_pomp_userdata_double (const char *name) {
  SEXP elt = getPairListElement(USERDATA,name);
  if (isNull(elt)) errorcall(R_NilValue,"no user-provided element named '%s' is found",name);
  if (!isReal(elt)) errorcall(R_NilValue,"user-provided element '%s' is not a numeric vector",name);
  return REAL(elt);
}

void unset_pomp_userdata (void) {
  USERDATA = R_NilValue;
}
