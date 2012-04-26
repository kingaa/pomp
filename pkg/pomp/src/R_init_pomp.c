#include "pomp.h"
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

void R_init_pomp (DllInfo *info) {
  R_RegisterCCallable("pomp","periodic_bspline_basis_eval",(DL_FUNC) &periodic_bspline_basis_eval);
  R_RegisterCCallable("pomp","dot_product",(DL_FUNC) &dot_product);
  R_RegisterCCallable("pomp","reulermultinom",(DL_FUNC) &reulermultinom);
  R_RegisterCCallable("pomp","deulermultinom",(DL_FUNC) &deulermultinom);
  R_RegisterCCallable("pomp","get_pomp_userdata",(DL_FUNC) &get_pomp_userdata);
}
