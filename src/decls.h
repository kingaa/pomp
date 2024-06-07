/* src/bspline.c */
extern SEXP bspline_basis(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP periodic_bspline_basis(SEXP, SEXP, SEXP, SEXP, SEXP);
extern void bspline_basis_eval_deriv(double, double *, int, int, int, double *);
extern void periodic_bspline_basis_eval_deriv(double, double, int, int, int, double *);
/* src/dinit.c */
extern SEXP do_dinit(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
/* src/distributions.c */
extern SEXP R_Euler_Multinom(SEXP, SEXP, SEXP, SEXP);
extern SEXP D_Euler_Multinom(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_GammaWN(SEXP, SEXP, SEXP);
extern SEXP R_BetaBinom(SEXP, SEXP, SEXP, SEXP);
extern SEXP D_BetaBinom(SEXP, SEXP, SEXP, SEXP, SEXP);
/* src/dmeasure.c */
extern SEXP do_dmeasure(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
/* src/dprior.c */
extern SEXP do_dprior(SEXP, SEXP, SEXP, SEXP);
/* src/dprocess.c */
extern SEXP do_dprocess(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
/* src/emeasure.c */
extern SEXP do_emeasure(SEXP, SEXP, SEXP, SEXP, SEXP);
/* src/euler.c */
extern SEXP euler_simulator(SEXP, SEXP, SEXP, SEXP, SEXP, double, rprocmode, SEXP, SEXP, SEXP, SEXP);
extern int num_euler_steps(double, double, double *);
extern int num_map_steps(double, double, double);
/* src/gompertz.c */
extern void _gompertz_normal_dmeasure(double *, double *, double *, double *, int, int *, int *, int *, int *, double *, double);
extern void _gompertz_normal_rmeasure(double *, double *, double *, int *, int *, int *, int *, double *, double);
extern void _gompertz_normal_emeasure(double *, double *, double *, int *, int *, int *, int *, double *, double);
extern void _gompertz_normal_vmeasure(double *, double *, double *, int *, int *, int *, int *, double *, double);
extern void _gompertz_step(double *, const double *, const int *, const int *, const int *, const double *, double, double);
extern void _gompertz_skeleton(double *, double *, const double *, const int *, const int *, const int *, const double *, double);
extern void _gompertz_to_trans(double *, const double *, const int *);
extern void _gompertz_from_trans(double *, const double *, const int *);
/* src/init.c */
extern void R_init_pomp(DllInfo *);
/* src/logmeanexp.c */
extern SEXP logmeanexp(const SEXP, const SEXP);
/* src/lookup_table.c */
extern SEXP get_covariate_names(SEXP);
extern lookup_table_t make_covariate_table(SEXP, int *);
extern SEXP lookup_in_table(SEXP, SEXP);
extern void table_lookup(lookup_table_t *, double, double *);
/* src/mif2.c */
extern SEXP randwalk_perturbation(SEXP, SEXP);
/* src/ou2.c */
extern void _ou2_step(double *, const double *, const int *, const int *, const int *, const double *, double, double);
extern void _ou2_pdf(double *, double *, double *, double, double, const double *, const int *, const int *, const int *, const double *);
extern void _ou2_skel(double *, double *, double *, int *, int *, int *, double *, double);
extern void _ou2_dmeasure(double *, double *, double *, double *, int, int *, int *, int *, int *, double *, double);
extern void _ou2_rmeasure(double *, double *, double *, int *, int *, int *, int *, double *, double);
extern void _ou2_emeasure(double *, double *, double *, int *, int *, int *, int *, double *, double);
extern void _ou2_vmeasure(double *, double *, double *, int *, int *, int *, int *, double *, double);
/* src/partrans.c */
extern SEXP do_partrans(SEXP, SEXP, SEXP, SEXP);
/* src/pfilter.c */
extern SEXP pfilter(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
/* src/pomp_fun.c */
extern SEXP pomp_fun_handler(SEXP, SEXP, pompfunmode *, SEXP, SEXP, SEXP, SEXP);
extern SEXP load_stack_incr(SEXP);
extern SEXP load_stack_decr(SEXP);
/* src/probe.c */
extern SEXP apply_probe_data(SEXP, SEXP);
extern SEXP apply_probe_sim(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
/* src/probe_acf.c */
extern SEXP probe_acf(SEXP, SEXP, SEXP);
extern SEXP probe_ccf(SEXP, SEXP, SEXP, SEXP);
/* src/probe_marginal.c */
extern SEXP probe_marginal_setup(SEXP, SEXP, SEXP);
extern SEXP probe_marginal_solve(SEXP, SEXP, SEXP);
/* src/probe_nlar.c */
extern SEXP probe_nlar(SEXP, SEXP, SEXP);
/* src/resample.c */
extern SEXP systematic_resampling(SEXP, SEXP);
extern void nosort_resamp(int, double *, int, int *, int);
/* src/rinit.c */
extern SEXP do_rinit(SEXP, SEXP, SEXP, SEXP, SEXP);
/* src/rmeasure.c */
extern SEXP do_rmeasure(SEXP, SEXP, SEXP, SEXP, SEXP);
/* src/rprior.c */
extern SEXP do_rprior(SEXP, SEXP, SEXP);
/* src/rprocess.c */
extern SEXP do_rprocess(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
/* src/simulate.c */
extern SEXP do_simulate(SEXP, SEXP, SEXP, SEXP, SEXP);
/* src/skeleton.c */
extern SEXP add_skel_args(SEXP, SEXP, SEXP, SEXP);
extern void eval_skeleton_R(double *, double *, double *, double *, SEXP, SEXP, SEXP, int, int, int, int, int, int, int, lookup_table_t *, double *);
extern void iterate_skeleton_R(double *, double, double, double *, double *, double *, SEXP, SEXP, SEXP, int, int, int, int, int, int, int, lookup_table_t *, int *, double *);
extern void eval_skeleton_native(double *, double *, double *, double *, int, int, int, int, int, int, int, int *, int *, int *, lookup_table_t *, pomp_skeleton *, SEXP, double *);
extern void iterate_skeleton_native(double *, double, double, double *, double *, double *, int, int, int, int, int, int, int, int *, int *, int *, lookup_table_t *, int *, pomp_skeleton *, SEXP, double *);
extern SEXP do_skeleton(SEXP, SEXP, SEXP, SEXP, SEXP);
/* src/sobolseq.c */
extern SEXP sobol_sequence(SEXP, SEXP);
/* src/ssa.c */
extern SEXP SSA_simulator(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
/* src/synth_lik.c */
extern SEXP synth_loglik(SEXP, SEXP);
/* src/trajectory.c */
extern SEXP iterate_map(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP pomp_desolve_setup(SEXP, SEXP, SEXP, SEXP);
extern void pomp_vf_eval(int *, double *, double *, double *, double *, int *);
extern SEXP pomp_desolve_takedown(void);
/* src/transformations.c */
extern SEXP LogitTransform(SEXP);
extern SEXP ExpitTransform(SEXP);
extern SEXP LogBarycentricTransform(SEXP);
extern SEXP InverseLogBarycentricTransform(SEXP);
/* src/userdata.c */
extern SEXP set_pomp_userdata(SEXP);
extern const SEXP get_userdata(const char *);
extern const int *get_userdata_int(const char *);
extern const double *get_userdata_double(const char *);
/* src/vmeasure.c */
extern SEXP do_vmeasure(SEXP, SEXP, SEXP, SEXP, SEXP);
/* src/wpfilter.c */
extern SEXP wpfilter(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
