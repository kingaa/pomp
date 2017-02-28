// -*- C++ -*-

#ifndef _POMP_INTERNAL_H_
#define _POMP_INTERNAL_H_

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

#include "pomp.h"

# define MATCHNAMES(X,N,W) (matchnames(GET_NAMES(X),(N),(W)))
# define MATCHROWNAMES(X,N,W) (matchnames(GET_ROWNAMES(GET_DIMNAMES(X)),(N),(W)))
# define MATCHCOLNAMES(X,N,W) (matchnames(GET_COLNAMES(GET_DIMNAMES(X)),(N),(W)))

// lookup-table structure, as used internally
typedef struct lookup_table {
  int length, width;
  int index;
  double *x;
  double *y;
} lookup_table;

// routine to compute number of discrete-time steps to take.
// used by plugins in 'euler.c' and map iterator in 'trajectory.c'
int num_map_steps (double t1, double t2, double dt);
int num_euler_steps (double t1, double t2, double *dt);

// simple linear interpolation of the lookup table (with derivative if desired)
void table_lookup (struct lookup_table *tab, double x, double *y);
struct lookup_table make_covariate_table (SEXP object, int *ncovars);
// setting dydt = 0 in the call to 'table_lookup' will bypass computation of the derivative
// void table_lookup (struct lookup_table *tab, double x, double *y, double *dydt);

// bspline.c
SEXP bspline_basis(SEXP x, SEXP degree, SEXP knots);
SEXP bspline_basis_function(SEXP x, SEXP i, SEXP degree, SEXP knots);
SEXP periodic_bspline_basis (SEXP x, SEXP nbasis, SEXP degree, SEXP period);

// dmeasure.c
SEXP do_dmeasure (SEXP object, SEXP y, SEXP x, SEXP times, SEXP params, SEXP log, SEXP gnsi);

// dprior.c
SEXP do_dprior (SEXP object, SEXP params, SEXP log, SEXP gnsi);

// dprocess.c
SEXP do_dprocess (SEXP object, SEXP x, SEXP times, SEXP params, SEXP log, SEXP gnsi);

// euler.c
SEXP euler_model_simulator (SEXP func, 
                            SEXP xstart, SEXP times, SEXP params, 
                            SEXP deltat, SEXP method, SEXP zeronames,
                            SEXP tcovar, SEXP covar, SEXP args, SEXP gnsi);
SEXP euler_model_density (SEXP func, 
			  SEXP x, SEXP times, SEXP params,
			  SEXP tcovar, SEXP covar, SEXP log, SEXP args, 
			  SEXP gnsi);

// eulermultinom.c
SEXP R_Euler_Multinom (SEXP n, SEXP size, SEXP rate, SEXP dt);
SEXP D_Euler_Multinom (SEXP x, SEXP size, SEXP rate, SEXP dt, SEXP log);
SEXP R_GammaWN (SEXP n, SEXP sigma, SEXP deltat);

// initstate.c
SEXP do_init_state (SEXP object, SEXP params, SEXP t0, SEXP nsim, SEXP gnsi);

// lookup_table.c
SEXP lookup_in_table (SEXP ttable, SEXP xtable, SEXP t);

// mif.c
SEXP mif_update (SEXP pfp, SEXP theta, SEXP gamma, SEXP varfactor, 
		 SEXP sigma, SEXP pars);

// mif2.c
SEXP randwalk_perturbation (SEXP params, SEXP rw_sd);

// ou2.c
void ou2_rmeasure (double *y, double *x, double *p, 
		   int *obsindex, int *stateindex, int *parindex, int *covindex,
		   int ncovar, double *covar, double t);
void ou2_dmeasure (double *lik, double *y, double *x, double *p, int give_log, 
		   int *obsindex, int *stateindex, int *parindex, int *covindex,
		   int covdim, double *covar, double t);
void ou2_step (double *x, const double *p,
	       const int *stateindex, const int *parindex, const int *covindex,
	       int ncovars, const double *covars,
	       double t, double dt);
void ou2_pdf (double *f, 
	      double *x1, double *x2, double t1, double t2, const double *p, 
	      const int *stateindex, const int *parindex, const int *covindex,
	      int ncovars, const double *covars);

// partrans.c
SEXP do_partrans (SEXP object, SEXP params, SEXP dir, SEXP gnsi);

// pfilter.c
SEXP pfilter_computations (SEXP x, SEXP params, SEXP Np,
			   SEXP rw_sd,
			   SEXP predmean, SEXP predvar,
			   SEXP filtmean, SEXP trackancestry, SEXP onepar,
			   SEXP weights, SEXP tol);

// pomp_fun.c
typedef enum {undef=-1,Rfun=0,native=1,regNative=2} pompfunmode;
SEXP pomp_fun_handler (SEXP pfun, SEXP gnsi, pompfunmode *mode);
SEXP load_stack_incr (SEXP pack);
SEXP load_stack_decr (SEXP pack);

// probe.c
SEXP apply_probe_data (SEXP object, SEXP probes);
SEXP apply_probe_sim (SEXP object, SEXP nsim, SEXP params, SEXP seed, SEXP probes, SEXP datval);

// probe_acf.c
SEXP probe_acf (SEXP x, SEXP lags, SEXP corr);
SEXP probe_ccf (SEXP x, SEXP y, SEXP lags, SEXP corr);

// probe_marginal.c
SEXP probe_marginal_setup (SEXP ref, SEXP order, SEXP diff);
SEXP probe_marginal_solve (SEXP x, SEXP setup, SEXP diff);

// probe_nlar.c
SEXP probe_nlar (SEXP x, SEXP lags, SEXP powers);

// resample.c
void nosort_resamp (int nw, double *w, int np, int *p, int offset);
SEXP systematic_resampling(SEXP weights);

// rmeasure.c
SEXP do_rmeasure (SEXP object, SEXP x, SEXP times, SEXP params, SEXP gnsi);

// rprior.c
SEXP do_rprior (SEXP object, SEXP params, SEXP gnsi);

// rprocess.c
SEXP do_rprocess (SEXP object, SEXP xstart, SEXP times, SEXP params, SEXP offset, SEXP gnsi);

// skeleton.c
SEXP do_skeleton (SEXP object, SEXP x, SEXP t, SEXP params, SEXP gnsi);
void eval_skeleton_native (double *f, double *time, double *x, double *p,
			   int nvars, int npars, int ncovars, int ntimes,
			   int nrepx, int nrepp, int nreps,
			   int *sidx, int *pidx, int *cidx,
			   lookup_table *covar_table,
			   pomp_skeleton *fun, SEXP args);
void eval_skeleton_R (double *f, double *time, double *x, double *p,
		      SEXP fcall, SEXP rho, SEXP Snames,
		      double *tp, double *xp, double *pp, double *cp,
		      int nvars, int npars, int ntimes,
		      int nrepx, int nrepp, int nreps, lookup_table *covar_table);

// simulate.c
SEXP simulation_computations (SEXP object, SEXP params, SEXP times, SEXP t0, 
			      SEXP nsim, SEXP obs, SEXP states, SEXP gnsi);

// sobolseq.c
SEXP sobol_sequence (SEXP dim, SEXP length);

// SSA_simulator.c
SEXP SSA_simulator (SEXP func, SEXP mflag, SEXP xstart, SEXP times, SEXP params,
		    SEXP e, SEXP vmatrix, SEXP dmatrix, SEXP deps, 
		    SEXP tcovar, SEXP covar,
		    SEXP zeronames, SEXP args, SEXP gnsi);

// synthlik.c
SEXP synth_loglik (SEXP ysim, SEXP ydat);

// trajectory.c
SEXP iterate_map (SEXP object, SEXP times, SEXP t0, SEXP x0, SEXP params, SEXP gnsi);
SEXP pomp_desolve_setup (SEXP object, SEXP x0, SEXP params, SEXP gnsi);
void pomp_desolve_takedown (void);

//userdata.c
void set_pomp_userdata (SEXP object);
void unset_pomp_userdata (void);

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
  PROTECT(index = match(x,names,0));
  idx = INTEGER(index);
  for (k = 0; k < n; k++) {
    if (idx[k]==0)
      errorcall(R_NilValue,"variable '%s' not found among the %s",
		CHARACTER_DATA(STRING_ELT(nm,k)),
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

#endif
