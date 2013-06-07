#include "pomp.h"
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

void R_init_pomp (DllInfo *info) {
  R_RegisterCCallable("pomp","periodic_bspline_basis_eval",(DL_FUNC) &periodic_bspline_basis_eval);
  R_RegisterCCallable("pomp","get_pomp_userdata",(DL_FUNC) &get_pomp_userdata);
  R_RegisterCCallable("pomp","get_pomp_userdata_int",(DL_FUNC) &get_pomp_userdata_int);
  R_RegisterCCallable("pomp","get_pomp_userdata_double",(DL_FUNC) &get_pomp_userdata_double);
}
