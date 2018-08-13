#include "pomp_internal.h"
#include <R_ext/Rdynload.h>

static const R_CallMethodDef callMethods[] = {
  {"bspline_basis", (DL_FUNC) &bspline_basis, 4},
  {"periodic_bspline_basis", (DL_FUNC) &periodic_bspline_basis, 5},
  {"systematic_resampling", (DL_FUNC) &systematic_resampling, 1},
  {"euler_model_simulator", (DL_FUNC) &euler_model_simulator, 11},
  {"onestep_density", (DL_FUNC) &onestep_density, 9},
  {"lookup_in_table", (DL_FUNC) &lookup_in_table, 2},
  {"load_stack_incr", (DL_FUNC) &load_stack_incr, 1},
  {"load_stack_decr", (DL_FUNC) &load_stack_decr, 1},
  {"SSA_simulator", (DL_FUNC) &SSA_simulator, 11},
  {"R_Euler_Multinom", (DL_FUNC) &R_Euler_Multinom, 4},
  {"D_Euler_Multinom", (DL_FUNC) &D_Euler_Multinom, 5},
  {"R_GammaWN", (DL_FUNC) &R_GammaWN, 3},
  {"pfilter_computations", (DL_FUNC) &pfilter_computations, 10},
  {"randwalk_perturbation", (DL_FUNC) &randwalk_perturbation, 2},
  {"simulation_computations", (DL_FUNC) &simulation_computations, 8},
  {"iterate_map", (DL_FUNC) &iterate_map, 6},
  {"pomp_desolve_setup", (DL_FUNC) &pomp_desolve_setup, 4},
  {"pomp_desolve_takedown", (DL_FUNC) &pomp_desolve_takedown, 0},
  {"apply_probe_data", (DL_FUNC) &apply_probe_data, 2},
  {"apply_probe_sim", (DL_FUNC) &apply_probe_sim, 7},
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
  {"do_rprior", (DL_FUNC) &do_rprior, 3},
  {"do_dprior", (DL_FUNC) &do_dprior, 4},
  {"do_skeleton", (DL_FUNC) &do_skeleton, 5},
  {"do_rinit", (DL_FUNC) &do_rinit, 5},
  {NULL, NULL, 0}
};

void R_init_pomp (DllInfo *info) {
  // C functions provided for users
  R_RegisterCCallable("pomp","periodic_bspline_basis_eval",(DL_FUNC) &periodic_bspline_basis_eval);
  R_RegisterCCallable("pomp","get_pomp_userdata",(DL_FUNC) &get_pomp_userdata);
  R_RegisterCCallable("pomp","get_pomp_userdata_int",(DL_FUNC) &get_pomp_userdata_int);
  R_RegisterCCallable("pomp","get_pomp_userdata_double",(DL_FUNC) &get_pomp_userdata_double);
  R_RegisterCCallable("pomp","pomp_fun_handler",(DL_FUNC) &pomp_fun_handler);
  R_RegisterCCallable("pomp","load_stack_incr",(DL_FUNC) &load_stack_incr);
  R_RegisterCCallable("pomp","load_stack_decr",(DL_FUNC) &load_stack_decr);
  R_RegisterCCallable("pomp","set_pomp_userdata",(DL_FUNC) &set_pomp_userdata);
  R_RegisterCCallable("pomp","unset_pomp_userdata",(DL_FUNC) &unset_pomp_userdata);
  R_RegisterCCallable("pomp","make_covariate_table",(DL_FUNC) &make_covariate_table);
  R_RegisterCCallable("pomp","table_lookup",(DL_FUNC) &table_lookup);
  R_RegisterCCallable("pomp","apply_probe_data",(DL_FUNC) &apply_probe_data);
  R_RegisterCCallable("pomp","apply_probe_sim",(DL_FUNC) &apply_probe_sim);
  R_RegisterCCallable("pomp","systematic_resampling",(DL_FUNC) &systematic_resampling);

  // Register routines
  R_registerRoutines(info,NULL,callMethods,NULL,NULL);
  R_useDynamicSymbols(info,TRUE);
}
