#include "internal.h"
#include <R_ext/Rdynload.h>

static const R_CallMethodDef callMethods[] = {
  {"set_userdata", (DL_FUNC) &set_pomp_userdata, 1},
  {"logmeanexp", (DL_FUNC) &logmeanexp, 2},
  {"bspline_basis", (DL_FUNC) &bspline_basis, 5},
  {"periodic_bspline_basis", (DL_FUNC) &periodic_bspline_basis, 5},
  {"systematic_resampling", (DL_FUNC) &systematic_resampling, 2},
  {"lookup_in_table", (DL_FUNC) &lookup_in_table, 2},
  {"load_stack_incr", (DL_FUNC) &load_stack_incr, 1},
  {"load_stack_decr", (DL_FUNC) &load_stack_decr, 1},
  {"R_Euler_Multinom", (DL_FUNC) &R_Euler_Multinom, 4},
  {"D_Euler_Multinom", (DL_FUNC) &D_Euler_Multinom, 5},
  {"E_Euler_Multinom", (DL_FUNC) &E_Euler_Multinom, 3},
  {"R_GammaWN", (DL_FUNC) &R_GammaWN, 3},
  {"R_BetaBinom", (DL_FUNC) &R_BetaBinom, 4},
  {"D_BetaBinom", (DL_FUNC) &D_BetaBinom, 5},
  {"pfilter_computations", (DL_FUNC) &pfilter, 10},
  {"wpfilter_comps", (DL_FUNC) &wpfilter, 7},
  {"randwalk_perturbation", (DL_FUNC) &randwalk_perturbation, 2},
  {"do_simulate", (DL_FUNC) &do_simulate, 5},
  {"iterate_map", (DL_FUNC) &iterate_map, 6},
  {"pomp_desolve_setup", (DL_FUNC) &pomp_desolve_setup, 4},
  {"pomp_desolve_takedown", (DL_FUNC) &pomp_desolve_takedown, 0},
  {"apply_probe_data", (DL_FUNC) &apply_probe_data, 2},
  {"apply_probe_sim", (DL_FUNC) &apply_probe_sim, 6},
  {"probe_marginal_setup", (DL_FUNC) &probe_marginal_setup, 3},
  {"probe_marginal_solve", (DL_FUNC) &probe_marginal_solve, 3},
  {"probe_acf", (DL_FUNC) &probe_acf, 3},
  {"probe_ccf", (DL_FUNC) &probe_ccf, 4},
  {"probe_nlar", (DL_FUNC) &probe_nlar, 3},
  {"synth_loglik", (DL_FUNC) &synth_loglik, 2},
  {"sobol_sequence", (DL_FUNC) &sobol_sequence, 2},
  {"do_partrans", (DL_FUNC) &do_partrans, 4},
  {"do_rprocess", (DL_FUNC) &do_rprocess, 6},
  {"do_dprocess", (DL_FUNC) &do_dprocess, 6},
  {"do_rmeasure", (DL_FUNC) &do_rmeasure, 5},
  {"do_dmeasure", (DL_FUNC) &do_dmeasure, 7},
  {"do_emeasure", (DL_FUNC) &do_emeasure, 5},
  {"do_vmeasure", (DL_FUNC) &do_vmeasure, 5},
  {"do_rprior", (DL_FUNC) &do_rprior, 3},
  {"do_dprior", (DL_FUNC) &do_dprior, 4},
  {"do_skeleton", (DL_FUNC) &do_skeleton, 5},
  {"do_rinit", (DL_FUNC) &do_rinit, 5},
  {"do_dinit", (DL_FUNC) &do_dinit, 6},
  {"LogitTransform", (DL_FUNC) &LogitTransform, 1},
  {"ExpitTransform", (DL_FUNC) &ExpitTransform, 1},
  {"LogBarycentricTransform", (DL_FUNC) &LogBarycentricTransform, 1},
  {"InverseLogBarycentricTransform", (DL_FUNC) &InverseLogBarycentricTransform, 1},
  {NULL, NULL, 0}
};

void R_init_pomp (DllInfo *info) {
  // C functions provided for users
  R_RegisterCCallable("pomp","bspline_basis_eval_deriv",(DL_FUNC) &bspline_basis_eval_deriv);
  R_RegisterCCallable("pomp","periodic_bspline_basis_eval_deriv",(DL_FUNC) &periodic_bspline_basis_eval_deriv);
  R_RegisterCCallable("pomp","get_userdata",(DL_FUNC) &get_userdata);
  R_RegisterCCallable("pomp","get_userdata_int",(DL_FUNC) &get_userdata_int);
  R_RegisterCCallable("pomp","get_userdata_double",(DL_FUNC) &get_userdata_double);
  R_RegisterCCallable("pomp","pomp_fun_handler",(DL_FUNC) &pomp_fun_handler);
  R_RegisterCCallable("pomp","load_stack_incr",(DL_FUNC) &load_stack_incr);
  R_RegisterCCallable("pomp","load_stack_decr",(DL_FUNC) &load_stack_decr);
  R_RegisterCCallable("pomp","make_covariate_table",(DL_FUNC) &make_covariate_table);
  R_RegisterCCallable("pomp","get_covariate_names",(DL_FUNC) &get_covariate_names);
  R_RegisterCCallable("pomp","table_lookup",(DL_FUNC) &table_lookup);
  R_RegisterCCallable("pomp","lookup_in_table",(DL_FUNC) &lookup_in_table);
  R_RegisterCCallable("pomp","apply_probe_data",(DL_FUNC) &apply_probe_data);
  R_RegisterCCallable("pomp","apply_probe_sim",(DL_FUNC) &apply_probe_sim);
  R_RegisterCCallable("pomp","systematic_resampling",(DL_FUNC) &systematic_resampling);
  R_RegisterCCallable("pomp","randwalk_perturbation", (DL_FUNC) &randwalk_perturbation);

  // Register routines
  R_registerRoutines(info,NULL,callMethods,NULL,NULL);
  R_useDynamicSymbols(info,TRUE);
  R_forceSymbols(info,FALSE);
}
