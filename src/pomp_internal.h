// -*- C++ -*-

#ifndef _POMP_INTERNAL_H_
#define _POMP_INTERNAL_H_

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

#include "pomp.h"

# define MATCHROWNAMES(X,N,W) (matchnames(GET_ROWNAMES(GET_DIMNAMES(X)),(N),(W)))
# define MATCHCOLNAMES(X,N,W) (matchnames(GET_COLNAMES(GET_DIMNAMES(X)),(N),(W)))

// lookup-table structure, as used internally
typedef struct lookup_table {
  int length, width;
  int index;
  double *x;
  double *y;
} lookup_table;

typedef enum {undef=-1,Rfun=0,native=1,regNative=2} pompfunmode;

static R_INLINE SEXP makearray (int rank, int *dim) {
  int *dimp, k;
  double *xp;
  SEXP dimx, x;
  PROTECT(dimx = NEW_INTEGER(rank));
  dimp = INTEGER(dimx);
  for (k = 0; k < rank; k++) dimp[k] = dim[k];
  PROTECT(x = allocArray(REALSXP,dimx));
  xp = REAL(x);
  for (k = 0; k < length(x); k++) xp[k] = NA_REAL;
  UNPROTECT(2);
  return x;
}

static R_INLINE SEXP matchnames (SEXP x, SEXP names, const char *where) {
  int n = length(names);
  int *idx, k;
  SEXP index, nm;
  PROTECT(nm = AS_CHARACTER(names));
  PROTECT(index = match(x,nm,0));
  idx = INTEGER(index);
  for (k = 0; k < n; k++) {
    if (idx[k]==0)
      errorcall(R_NilValue,"variable '%s' not found among the %s",
		CHAR(STRING_ELT(nm,k)),
		where);
    idx[k] -= 1;
  }
  UNPROTECT(2);
  return index;
}

static R_INLINE SEXP name_index (SEXP names, SEXP object, const char *slot, const char *humanreadable) {
  SEXP slotnames, index;
  PROTECT(slotnames = GET_SLOT(object,install(slot)));
  if (LENGTH(slotnames) > 0) {
    PROTECT(index = matchnames(names,slotnames,humanreadable));
  } else {
    PROTECT(index = NEW_INTEGER(0));
  }
  UNPROTECT(2);
  return index;
}

static R_INLINE void setrownames (SEXP x, SEXP names, int n) {
  SEXP dimnms, nm;
  PROTECT(nm = AS_CHARACTER(names));
  PROTECT(dimnms = allocVector(VECSXP,n));
  SET_ELEMENT(dimnms,0,nm);	// set row names
  SET_DIMNAMES(x,dimnms);
  UNPROTECT(2);
}

static R_INLINE void fixdimnames (SEXP x, const char **names, int n) {
  int nprotect = 0;
  int i;
  SEXP dimnames, nm;
  PROTECT(dimnames = GET_DIMNAMES(x)); nprotect++;
  if (isNull(dimnames)) {
    PROTECT(dimnames = allocVector(VECSXP,n)); nprotect++;
  }
  PROTECT(nm = allocVector(VECSXP,n)); nprotect++;
  for (i = 0; i < n; i++)
    SET_ELEMENT(nm,i,mkChar(names[i]));
  SET_NAMES(dimnames,nm);
  SET_DIMNAMES(x,dimnames);
  UNPROTECT(nprotect);
}

static R_INLINE SEXP as_matrix (SEXP x) {
  int nprotect = 0;
  SEXP dim, names;
  int *xdim, nrow, ncol;
  PROTECT(dim = GET_DIM(x)); nprotect++;
  if (isNull(dim)) {
    PROTECT(x = duplicate(x)); nprotect++;
    PROTECT(names = GET_NAMES(x)); nprotect++;
    dim = NEW_INTEGER(2);
    xdim = INTEGER(dim);
    xdim[0] = LENGTH(x); xdim[1] = 1;
    SET_DIM(x,dim);
    SET_NAMES(x,R_NilValue);
    setrownames(x,names,2);
  } else if (LENGTH(dim) == 1) {
    PROTECT(x = duplicate(x)); nprotect++;
    PROTECT(names = GET_ROWNAMES(GET_DIMNAMES(x))); nprotect++;
    dim = NEW_INTEGER(2);
    xdim = INTEGER(dim);
    xdim[0] = LENGTH(x); xdim[1] = 1;
    SET_DIM(x,dim);
    SET_NAMES(x,R_NilValue);
    setrownames(x,names,2);
  } else if (LENGTH(dim) > 2) {
    PROTECT(x = duplicate(x)); nprotect++;
    PROTECT(names = GET_ROWNAMES(GET_DIMNAMES(x))); nprotect++;
    nrow = INTEGER(dim)[0];
    ncol = LENGTH(x)/nrow;
    dim = NEW_INTEGER(2);
    xdim = INTEGER(dim);
    xdim[0] = nrow; xdim[1] = ncol;
    SET_DIM(x,dim);
    SET_NAMES(x,R_NilValue);
    setrownames(x,names,2);
  }
  UNPROTECT(nprotect);
  return x;
}

static R_INLINE SEXP as_state_array (SEXP x) {
  int nprotect = 0;
  SEXP dim, names;
  int *xdim, nrow, ncol;
  PROTECT(dim = GET_DIM(x)); nprotect++;
  if (isNull(dim)) {
    PROTECT(x = duplicate(x)); nprotect++;
    PROTECT(names = GET_NAMES(x)); nprotect++;
    dim = NEW_INTEGER(3);
    xdim = INTEGER(dim);
    xdim[0] = LENGTH(x); xdim[1] = 1; xdim[2] = 1;
    SET_DIM(x,dim);
    SET_NAMES(x,R_NilValue);
    setrownames(x,names,3);
  } else if (LENGTH(dim) == 1) {
    PROTECT(x = duplicate(x)); nprotect++;
    PROTECT(names = GET_ROWNAMES(GET_DIMNAMES(x))); nprotect++;
    dim = NEW_INTEGER(3);
    xdim = INTEGER(dim);
    xdim[0] = LENGTH(x); xdim[1] = 1; xdim[2] = 1;
    SET_DIM(x,dim);
    SET_NAMES(x,R_NilValue);
    setrownames(x,names,3);
  } else if (LENGTH(dim) == 2) {
    PROTECT(x = duplicate(x)); nprotect++;
    PROTECT(names = GET_ROWNAMES(GET_DIMNAMES(x))); nprotect++;
    xdim = INTEGER(dim);
    nrow = xdim[0]; ncol = xdim[1];
    dim = NEW_INTEGER(3);
    xdim = INTEGER(dim);
    xdim[0] = nrow; xdim[1] = 1; xdim[2] = ncol;
    SET_DIM(x,dim);
    SET_NAMES(x,R_NilValue);
    setrownames(x,names,3);
  } else if (LENGTH(dim) > 3) {
    PROTECT(x = duplicate(x)); nprotect++;
    PROTECT(names = GET_ROWNAMES(GET_DIMNAMES(x))); nprotect++;
    xdim = INTEGER(dim);
    nrow = xdim[0]; ncol = xdim[1];
    dim = NEW_INTEGER(3);
    xdim = INTEGER(dim);
    xdim[0] = nrow; xdim[1] = ncol; xdim[2] = LENGTH(x)/nrow/ncol;
    SET_DIM(x,dim);
    SET_NAMES(x,R_NilValue);
    setrownames(x,names,3);
  }
  UNPROTECT(nprotect);
  return x;
}

static R_INLINE SEXP getListElement (SEXP list, const char *str)
{
  SEXP elmt = R_NilValue;
  SEXP names = getAttrib(list,R_NamesSymbol);
  for (R_len_t i = 0; i < length(list); i++)
    if(strcmp(CHAR(STRING_ELT(names,i)),str) == 0) {
      elmt = VECTOR_ELT(list,i);
      break;
    }
  return elmt;
}

static R_INLINE SEXP getPairListElement (SEXP list, const char *name)
{
  const char *tag;
  while (list != R_NilValue) {
    tag = CHAR(PRINTNAME(TAG(list)));
    if (strcmp(tag,name)==0) break;
    list = CDR(list);
  }
  return CAR(list);
}

#ifdef __cplusplus

template <class Scalar>
class view {
private:
  Scalar *data;
  int dim[2];
public:
  view (Scalar *x) {
    data = x;
    dim[0] = 0;
    dim[1] = 0;
  };
  view (Scalar *x, int d1) {
    data = x;
    dim[0] = d1;
    dim[1] = 0;
  };
  view (Scalar *x, int d1, int d2) {
    data = x;
    dim[0] = d1;
    dim[1] = d2;
  };
  ~view (void) {};
  inline Scalar& operator () (int d1) {
    return(data[d1]);
  };
  inline Scalar& operator () (int d1, int d2) {
    return(data[d1 + dim[0] * d2]);
  };
  inline Scalar& operator () (int d1, int d2, int d3) {
    return(data[d1 + dim[0] * (d2 + dim[1] * d3)]);
  };
};

#endif

// PROTOTYPES
// Make with, e.g.,
// cproto -I `R RHOME`/include -e *.c | perl -pe 's/\/\*(.+?)\*\//\n\/\/$1/g'

// blowfly.c
extern void _blowfly_dmeasure(double *lik, double *y, double *x, double *p, int give_log, int *obsindex, int *stateindex, int *parindex, int *covindex, int ncovars, double *covars, double t);
extern void _blowfly_rmeasure(double *y, double *x, double *p, int *obsindex, int *stateindex, int *parindex, int *covindex, int ncovars, double *covars, double t);
extern void _blowfly_simulator_one(double *x, const double *p, const int *stateindex, const int *parindex, const int *covindex, int covdim, const double *covar, double t, double Rf_dt);
extern void _blowfly_simulator_two(double *x, const double *p, const int *stateindex, const int *parindex, const int *covindex, int covdim, const double *covar, double t, double Rf_dt);

// bspline.c
extern SEXP bspline_basis(SEXP x, SEXP nbasis, SEXP degree, SEXP deriv);
extern SEXP periodic_bspline_basis(SEXP x, SEXP nbasis, SEXP degree, SEXP period, SEXP deriv);
extern void periodic_bspline_basis_eval(double x, double period, int degree, int nbasis, double *y);
extern void periodic_bspline_basis_eval_deriv(double x, double period, int degree, int nbasis, int deriv, double *y);

// dmeasure.c
extern SEXP do_dmeasure(SEXP object, SEXP y, SEXP x, SEXP times, SEXP params, SEXP log, SEXP gnsi);

// dprior.c
extern void _pomp_default_dprior(double *lik, double *p, int give_log, int *parindex);
extern SEXP do_dprior(SEXP object, SEXP params, SEXP log, SEXP gnsi);

// dprocess.c
extern SEXP euler_model_density(SEXP func, SEXP x, SEXP times, SEXP params, SEXP tcovar, SEXP covar, SEXP log, SEXP args, SEXP gnsi);
extern SEXP do_dprocess(SEXP object, SEXP x, SEXP times, SEXP params, SEXP log, SEXP gnsi);

// euler.c
extern SEXP euler_model_simulator(SEXP func, SEXP xstart, SEXP times, SEXP params, double deltat, int method, SEXP zeronames, SEXP tcovar, SEXP covar, SEXP args, SEXP gnsi);
extern int num_euler_steps(double t1, double t2, double *Rf_dt);
extern int num_map_steps(double t1, double t2, double Rf_dt);

// eulermultinom.c
extern SEXP R_Euler_Multinom(SEXP n, SEXP size, SEXP rate, SEXP Rf_dt);
extern SEXP D_Euler_Multinom(SEXP x, SEXP size, SEXP rate, SEXP Rf_dt, SEXP log);
extern SEXP R_GammaWN(SEXP n, SEXP sigma, SEXP deltat);

// gompertz.c
extern void _gompertz_normal_dmeasure(double *lik, double *y, double *x, double *p, int give_log, int *obsindex, int *stateindex, int *parindex, int *covindex, int ncovars, double *covars, double t);
extern void _gompertz_normal_rmeasure(double *y, double *x, double *p, int *obsindex, int *stateindex, int *parindex, int *covindex, int ncovars, double *covars, double t);
extern void _gompertz_simulator(double *x, const double *p, const int *stateindex, const int *parindex, const int *covindex, int covdim, const double *covar, double t, double Rf_dt);
extern void _gompertz_skeleton(double *f, double *x, const double *p, const int *stateindex, const int *parindex, const int *covindex, int covdim, const double *covar, double t);

// init.c
extern void R_init_pomp(DllInfo *info);

// initstate.c
extern SEXP do_init_state(SEXP object, SEXP params, SEXP t0, SEXP nsim, SEXP gnsi);

// lookup_table.c
extern struct lookup_table make_covariate_table(SEXP object, int *ncovar);
extern SEXP lookup_in_table(SEXP ttable, SEXP xtable, SEXP t);
extern void table_lookup(struct lookup_table *tab, double x, double *y);

// mif2.c
extern SEXP randwalk_perturbation(SEXP params, SEXP rw_sd);

// mif.c
extern SEXP mif_update(SEXP pfp, SEXP theta, SEXP gamma, SEXP varfactor, SEXP sigma, SEXP pars);

// ou2.c
extern void _ou2_step(double *x, const double *p, const int *stateindex, const int *parindex, const int *covindex, int ncovars, const double *covars, double t, double Rf_dt);
extern void _ou2_pdf(double *f, double *x, double *z, double t1, double t2, const double *p, const int *stateindex, const int *parindex, const int *covindex, int ncovars, const double *covars);
extern void _ou2_skel(double *f, double *x, double *p, int *stateindex, int *parindex, int *covindex, int ncovars, double *covars, double t);
extern void _ou2_dmeasure(double *lik, double *y, double *x, double *p, int give_log, int *obsindex, int *stateindex, int *parindex, int *covindex, int covdim, double *covar, double t);
extern void _ou2_rmeasure(double *y, double *x, double *p, int *obsindex, int *stateindex, int *parindex, int *covindex, int ncovar, double *covar, double t);

// partrans.c
extern SEXP do_partrans(SEXP object, SEXP params, SEXP dir, SEXP gnsi);

// pfilter.c
extern SEXP pfilter_computations(SEXP x, SEXP params, SEXP Np, SEXP rw_sd, SEXP predmean, SEXP predvar, SEXP filtmean, SEXP trackancestry, SEXP onepar, SEXP weights, SEXP tol);

// pomp_fun.c
extern SEXP pomp_fun_handler(SEXP pfun, SEXP gnsi, pompfunmode *mode);
extern SEXP load_stack_incr(SEXP pack);
extern SEXP load_stack_decr(SEXP pack);

// probe_acf.c
extern SEXP probe_acf(SEXP x, SEXP lags, SEXP corr);
extern SEXP probe_ccf(SEXP x, SEXP y, SEXP lags, SEXP corr);

// probe.c
extern SEXP apply_probe_data(SEXP object, SEXP probes);
extern SEXP apply_probe_sim(SEXP object, SEXP nsim, SEXP params, SEXP seed, SEXP probes, SEXP datval);

// probe_marginal.c
extern SEXP probe_marginal_setup(SEXP ref, SEXP order, SEXP diff);
extern SEXP probe_marginal_solve(SEXP x, SEXP setup, SEXP diff);

// probe_nlar.c
extern SEXP probe_nlar(SEXP x, SEXP lags, SEXP powers);

// resample.c
extern void nosort_resamp(int nw, double *w, int np, int *p, int offset);
extern SEXP systematic_resampling(SEXP weights);

// ricker.c
extern void _ricker_poisson_dmeasure(double *lik, double *y, double *x, double *p, int give_log, int *obsindex, int *stateindex, int *parindex, int *covindex, int ncovars, double *covars, double t);
extern void _ricker_poisson_rmeasure(double *y, double *x, double *p, int *obsindex, int *stateindex, int *parindex, int *covindex, int ncovars, double *covars, double t);
extern void _ricker_simulator(double *x, const double *p, const int *stateindex, const int *parindex, const int *covindex, int covdim, const double *covar, double t, double Rf_dt);
extern void _ricker_skeleton(double *f, double *x, const double *p, const int *stateindex, const int *parindex, const int *covindex, int covdim, const double *covar, double t);

// rmeasure.c
extern SEXP do_rmeasure(SEXP object, SEXP x, SEXP times, SEXP params, SEXP gnsi);

// rprior.c
extern SEXP do_rprior(SEXP object, SEXP params, SEXP gnsi);

// rprocess.c
extern SEXP do_rprocess(SEXP object, SEXP xstart, SEXP times, SEXP params, SEXP offset, SEXP gnsi);

// simulate.c
extern SEXP simulation_computations(SEXP object, SEXP params, SEXP times, SEXP t0, SEXP nsim, SEXP obs, SEXP states, SEXP gnsi);

// sir.c
extern void _sir_par_untrans(double *Rf_pt, double *p, int *parindex);
extern void _sir_par_trans(double *Rf_pt, double *p, int *parindex);
extern void _sir_init(double *x, const double *p, double t, const int *stateindex, const int *parindex, const int *covindex, const double *covars);
extern void _sir_binom_dmeasure(double *lik, double *y, double *x, double *p, int give_log, int *obsindex, int *stateindex, int *parindex, int *covindex, int ncovars, double *covars, double t);
extern void _sir_binom_rmeasure(double *y, double *x, double *p, int *obsindex, int *stateindex, int *parindex, int *covindex, int ncovars, double *covars, double t);
extern void _sir_euler_simulator(double *x, const double *p, const int *stateindex, const int *parindex, const int *covindex, int covdim, const double *covar, double t, double Rf_dt);
extern void _sir_ODE(double *f, double *x, const double *p, const int *stateindex, const int *parindex, const int *covindex, int covdim, const double *covar, double t);
extern double _sir_rates(int j, double t, double *x, double *p, int *stateindex, int *parindex, int *covindex, int ncovar, double *covar);

// skeleton.c
extern void eval_skeleton_native(double *f, double *time, double *x, double *p, int nvars, int npars, int ncovars, int ntimes, int nrepx, int nrepp, int nreps, int *sidx, int *pidx, int *cidx, lookup_table *covar_table, pomp_skeleton *fun, SEXP args);
extern void eval_skeleton_R(double *f, double *time, double *x, double *p, SEXP fcall, SEXP rho, SEXP Snames, double *tp, double *xp, double *pp, double *cp, int nvars, int npars, int ntimes, int nrepx, int nrepp, int nreps, lookup_table *covar_table);
extern SEXP do_skeleton(SEXP object, SEXP x, SEXP t, SEXP params, SEXP gnsi);

// sobolseq.c
extern SEXP sobol_sequence(SEXP dim, SEXP length);

// ssa.c
extern SEXP SSA_simulator(SEXP func, SEXP xstart, SEXP times, SEXP params, SEXP vmatrix, SEXP tcovar, SEXP covar, SEXP zeronames, SEXP hmax, SEXP args, SEXP gnsi);

// synth_lik.c
extern SEXP synth_loglik(SEXP ysim, SEXP ydat);

// trajectory.c
extern SEXP iterate_map(SEXP object, SEXP times, SEXP t0, SEXP x0, SEXP params, SEXP gnsi);
extern SEXP pomp_desolve_setup(SEXP object, SEXP x0, SEXP params, SEXP gnsi);
extern void pomp_vf_eval(int *neq, double *t, double *y, double *ydot, double *yout, int *ip);
extern void pomp_desolve_takedown(void);

// userdata.c
extern void set_pomp_userdata(SEXP userdata);
extern const SEXP get_pomp_userdata(const char *name);
extern const int *get_pomp_userdata_int(const char *name);
extern const double *get_pomp_userdata_double(const char *name);
extern void unset_pomp_userdata(void);

#endif
