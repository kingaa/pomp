/* src/bspline.c */
extern void bspline_eval(double *y, const double *x, int nx, int i, int degree, int deriv, const double *knots);
extern SEXP bspline_basis(SEXP range, SEXP x, SEXP nbasis, SEXP degree, SEXP deriv);
extern SEXP periodic_bspline_basis(SEXP x, SEXP nbasis, SEXP degree, SEXP period, SEXP deriv);
extern void periodic_bspline_basis_eval(double x, double period, int degree, int nbasis, double *y);
extern void periodic_bspline_basis_eval_deriv(double x, double period, int degree, int nbasis, int deriv, double *y);
/* src/dinit.c */
extern SEXP do_dinit(SEXP object, SEXP t0, SEXP x, SEXP params, SEXP log, SEXP gnsi);
/* src/distributions.c */
extern SEXP R_Euler_Multinom(SEXP n, SEXP size, SEXP rate, SEXP deltat);
extern SEXP D_Euler_Multinom(SEXP x, SEXP size, SEXP rate, SEXP deltat, SEXP log);
extern SEXP R_GammaWN(SEXP n, SEXP sigma, SEXP deltat);
extern SEXP R_BetaBinom(SEXP n, SEXP size, SEXP prob, SEXP theta);
extern SEXP D_BetaBinom(SEXP x, SEXP size, SEXP prob, SEXP theta, SEXP log);
/* src/dmeasure.c */
extern SEXP do_dmeasure(SEXP object, SEXP y, SEXP x, SEXP times, SEXP params, SEXP log, SEXP gnsi);
/* src/dprior.c */
extern SEXP do_dprior(SEXP object, SEXP params, SEXP log, SEXP gnsi);
/* src/dprocess.c */
extern SEXP do_dprocess(SEXP object, SEXP x, SEXP times, SEXP params, SEXP log, SEXP gnsi);
/* src/emeasure.c */
extern SEXP do_emeasure(SEXP object, SEXP x, SEXP times, SEXP params, SEXP gnsi);
/* src/euler.c */
extern SEXP euler_model_simulator(SEXP func, SEXP xstart, SEXP tstart, SEXP times, SEXP params, double deltat, rprocmode method, SEXP accumvars, SEXP covar, SEXP args, SEXP gnsi);
extern int num_euler_steps(double t1, double t2, double *deltat);
extern int num_map_steps(double t1, double t2, double deltat);
/* src/gompertz.c */
extern void _gompertz_normal_dmeasure(double *lik, double *y, double *x, double *p, int give_log, int *obsindex, int *stateindex, int *parindex, int *covindex, double *covars, double t);
extern void _gompertz_normal_rmeasure(double *y, double *x, double *p, int *obsindex, int *stateindex, int *parindex, int *covindex, double *covars, double t);
extern void _gompertz_normal_emeasure(double *f, double *x, double *p, int *obsindex, int *stateindex, int *parindex, int *covindex, double *covars, double t);
extern void _gompertz_normal_vmeasure(double *f, double *x, double *p, int *vmatindex, int *stateindex, int *parindex, int *covindex, double *covars, double t);
extern void _gompertz_step(double *x, const double *p, const int *stateindex, const int *parindex, const int *covindex, const double *covar, double t, double deltat);
extern void _gompertz_skeleton(double *f, double *x, const double *p, const int *stateindex, const int *parindex, const int *covindex, const double *covar, double t);
extern void _gompertz_to_trans(double *__pt, const double *__p, const int *__parindex);
extern void _gompertz_from_trans(double *__p, const double *__pt, const int *__parindex);
/* src/init.c */
extern void R_init_pomp(DllInfo *info);
/* src/logmeanexp.c */
extern SEXP logmeanexp(const SEXP X, const SEXP Drop);
/* src/lookup_table.c */
extern SEXP get_covariate_names(SEXP object);
extern lookup_table_t make_covariate_table(SEXP object, int *ncovar);
extern SEXP lookup_in_table(SEXP covar, SEXP t);
extern void table_lookup(lookup_table_t *tab, double x, double *y);
/* src/mif2.c */
extern SEXP randwalk_perturbation(SEXP params, SEXP rw_sd);
/* src/ou2.c */
extern void _ou2_step(double *x, const double *p, const int *stateindex, const int *parindex, const int *covindex, const double *covars, double t, double deltat);
extern void _ou2_pdf(double *f, double *x, double *z, double t1, double t2, const double *p, const int *stateindex, const int *parindex, const int *covindex, const double *covars);
extern void _ou2_skel(double *f, double *x, double *p, int *stateindex, int *parindex, int *covindex, double *covars, double t);
extern void _ou2_dmeasure(double *lik, double *y, double *x, double *p, int give_log, int *obsindex, int *stateindex, int *parindex, int *covindex, double *covar, double t);
extern void _ou2_rmeasure(double *y, double *x, double *p, int *obsindex, int *stateindex, int *parindex, int *covindex, double *covar, double t);
extern void _ou2_emeasure(double *y, double *x, double *p, int *obsindex, int *stateindex, int *parindex, int *covindex, double *covar, double t);
extern void _ou2_vmeasure(double *f, double *x, double *p, int *vmatindex, int *stateindex, int *parindex, int *covindex, double *covar, double t);
/* src/partrans.c */
extern SEXP do_partrans(SEXP object, SEXP params, SEXP dir, SEXP gnsi);
/* src/pfilter.c */
extern SEXP pfilter(SEXP x, SEXP params, SEXP Np, SEXP predmean, SEXP predvar, SEXP filtmean, SEXP trackancestry, SEXP doparRS, SEXP weights, SEXP wave);
/* src/pomp_fun.c */
extern SEXP pomp_fun_handler(SEXP pfun, SEXP gnsi, pompfunmode *mode, SEXP S, SEXP P, SEXP O, SEXP C);
extern SEXP load_stack_incr(SEXP pack);
extern SEXP load_stack_decr(SEXP pack);
/* src/probe.c */
extern SEXP apply_probe_data(SEXP object, SEXP probes);
extern SEXP apply_probe_sim(SEXP object, SEXP nsim, SEXP params, SEXP probes, SEXP datval, SEXP gnsi);
/* src/probe_acf.c */
extern SEXP probe_acf(SEXP x, SEXP lags, SEXP corr);
extern SEXP probe_ccf(SEXP x, SEXP y, SEXP lags, SEXP corr);
/* src/probe_marginal.c */
extern SEXP probe_marginal_setup(SEXP ref, SEXP order, SEXP diff);
extern SEXP probe_marginal_solve(SEXP x, SEXP setup, SEXP diff);
/* src/probe_nlar.c */
extern SEXP probe_nlar(SEXP x, SEXP lags, SEXP powers);
/* src/resample.c */
extern SEXP systematic_resampling(SEXP weights, SEXP np);
extern void nosort_resamp(int nw, double *w, int np, int *p, int offset);
/* src/rinit.c */
extern SEXP do_rinit(SEXP object, SEXP params, SEXP t0, SEXP nsim, SEXP gnsi);
/* src/rmeasure.c */
extern SEXP do_rmeasure(SEXP object, SEXP x, SEXP times, SEXP params, SEXP gnsi);
/* src/rprior.c */
extern SEXP do_rprior(SEXP object, SEXP params, SEXP gnsi);
/* src/rprocess.c */
extern SEXP do_rprocess(SEXP object, SEXP xstart, SEXP tstart, SEXP times, SEXP params, SEXP gnsi);
/* src/simulate.c */
extern SEXP do_simulate(SEXP object, SEXP params, SEXP nsim, SEXP rettype, SEXP gnsi);
/* src/skeleton.c */
extern SEXP add_skel_args(SEXP args, SEXP Snames, SEXP Pnames, SEXP Cnames);
extern void eval_skeleton_R(double *f, double *time, double *x, double *p, SEXP fn, SEXP args, SEXP Snames, int nvars, int npars, int ncovars, int ntimes, int nrepx, int nrepp, int nreps, lookup_table_t *covar_table, double *cov);
extern void iterate_skeleton_R(double *X, double t, double deltat, double *time, double *x, double *p, SEXP fn, SEXP args, SEXP Snames, int nvars, int npars, int ncovars, int ntimes, int nrepp, int nreps, int nzeros, lookup_table_t *covar_table, int *zeroindex, double *cov);
extern void eval_skeleton_native(double *f, double *time, double *x, double *p, int nvars, int npars, int ncovars, int ntimes, int nrepx, int nrepp, int nreps, int *sidx, int *pidx, int *cidx, lookup_table_t *covar_table, pomp_skeleton *fun, SEXP args, double *cov);
extern void iterate_skeleton_native(double *X, double t, double deltat, double *time, double *x, double *p, int nvars, int npars, int ncovars, int ntimes, int nrepp, int nreps, int nzeros, int *sidx, int *pidx, int *cidx, lookup_table_t *covar_table, int *zeroindex, pomp_skeleton *fun, SEXP args, double *cov);
extern SEXP do_skeleton(SEXP object, SEXP x, SEXP t, SEXP params, SEXP gnsi);
/* src/sobolseq.c */
extern SEXP sobol_sequence(SEXP dim, SEXP length);
/* src/ssa.c */
extern SEXP SSA_simulator(SEXP func, SEXP xstart, SEXP tstart, SEXP times, SEXP params, SEXP vmatrix, SEXP covar, SEXP accumvars, SEXP hmax, SEXP args, SEXP gnsi);
/* src/synth_lik.c */
extern SEXP synth_loglik(SEXP ysim, SEXP ydat);
/* src/trajectory.c */
extern SEXP iterate_map(SEXP object, SEXP times, SEXP t0, SEXP x0, SEXP params, SEXP gnsi);
extern SEXP pomp_desolve_setup(SEXP object, SEXP x0, SEXP params, SEXP gnsi);
extern void pomp_vf_eval(int *neq, double *t, double *y, double *ydot, double *yout, int *ip);
extern void pomp_desolve_takedown(void);
/* src/transformations.c */
extern SEXP LogitTransform(SEXP P);
extern SEXP ExpitTransform(SEXP X);
extern SEXP LogBarycentricTransform(SEXP X);
extern SEXP InverseLogBarycentricTransform(SEXP Y);
/* src/userdata.c */
extern void set_pomp_userdata(SEXP userdata);
extern const SEXP get_userdata(const char *name);
extern const int *get_userdata_int(const char *name);
extern const double *get_userdata_double(const char *name);
extern void unset_pomp_userdata(void);
/* src/vmeasure.c */
extern SEXP do_vmeasure(SEXP object, SEXP x, SEXP times, SEXP params, SEXP gnsi);
/* src/wpfilter.c */
extern SEXP wpfilter(SEXP X, SEXP Params, SEXP Weights, SEXP W, SEXP Trigger, SEXP Target, SEXP Np);
