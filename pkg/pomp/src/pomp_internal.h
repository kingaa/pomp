// -*- C++ -*-

#ifndef _POMP_INTERNAL_H_
#define _POMP_INTERNAL_H_

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

#include "pomp.h"

# define MATCHROWNAMES(X,N) (matchnames(GET_ROWNAMES(GET_DIMNAMES(X)),(N)))
# define MATCHCOLNAMES(X,N) (matchnames(GET_COLNAMES(GET_DIMNAMES(X)),(N)))
# define MATCH_CHAR_TO_ROWNAMES(X,N,A) (match_char_to_names(GET_ROWNAMES(GET_DIMNAMES(X)),(N),(A)))

// lookup-table structure, as used internally
struct lookup_table {
  int length, width;
  int index;
  double *x;
  double *y;
};

// routine to compute number of discrete-time steps to take.
// used by plugins in 'euler.c' and map iterator in 'trajectory.c'
int num_map_steps (double t1, double t2, double dt);

// simple linear interpolation of the lookup table (with derivative if desired)
// setting dydt = 0 in the call to 'table_lookup' will bypass computation of the derivative
void table_lookup (struct lookup_table *tab, double x, double *y, double *dydt);

// bspline.c
SEXP bspline_basis(SEXP x, SEXP degree, SEXP knots);
SEXP bspline_basis_function(SEXP x, SEXP i, SEXP degree, SEXP knots);

// dsobol.c
SEXP sobol_sequence(SEXP dim);

// pomp_fun.c
SEXP pomp_fun_handler (SEXP pfun, int *use_native);
SEXP get_pomp_fun (SEXP pfun);

// lookup_table.c
SEXP lookup_in_table (SEXP ttable, SEXP xtable, SEXP t);

// resample.c
SEXP systematic_resampling(SEXP weights);

// initstate.c
SEXP do_init_state (SEXP object, SEXP params, SEXP t0);

// rprocess.c
SEXP do_rprocess (SEXP object, SEXP xstart, SEXP times, SEXP params, SEXP offset);

// rmeasure.c
SEXP do_rmeasure (SEXP object, SEXP x, SEXP times, SEXP params, SEXP fun);

// skeleton.c
SEXP do_skeleton (SEXP object, SEXP x, SEXP t, SEXP params, SEXP fun);

//userdata.c
void set_pomp_userdata (SEXP object);
void unset_pomp_userdata (void);

static R_INLINE SEXP makearray (int rank, int *dim) {
  int nprotect = 0;
  int *dimp, k;
  double *xp;
  SEXP dimx, x;
  PROTECT(dimx = NEW_INTEGER(rank)); nprotect++;
  dimp = INTEGER(dimx); 
  for (k = 0; k < rank; k++) dimp[k] = dim[k];
  PROTECT(x = allocArray(REALSXP,dimx)); nprotect++;
  xp = REAL(x);
  for (k = 0; k < length(x); k++) xp[k] = NA_REAL;
  UNPROTECT(nprotect);
  return x;
}

static R_INLINE SEXP matchnames (SEXP x, SEXP names) {
  int nprotect = 0;
  int n = length(names);
  int *idx, k;
  SEXP index, nm;
  PROTECT(nm = AS_CHARACTER(names)); nprotect++;
  PROTECT(index = match(x,names,0)); nprotect++;
  idx = INTEGER(index);
  for (k = 0; k < n; k++) {
    if (idx[k]==0) error("variable '%s' not found",CHARACTER_DATA(STRING_ELT(nm,k)));
    idx[k] -= 1;
  }
  UNPROTECT(nprotect);
  return index;
}

static R_INLINE SEXP match_char_to_names (SEXP x, int n, char **names) {
  int nprotect = 0;
  int *idx, k;
  SEXP index, nm;
  PROTECT(nm = NEW_CHARACTER(n)); nprotect++;
  for (k = 0; k < n; k++) {
    SET_STRING_ELT(nm,k,mkChar(names[k]));
  }
  PROTECT(index = match(x,nm,0)); nprotect++;
  idx = INTEGER(index);
  for (k = 0; k < n; k++) {
    if (idx[k]==0) {
      UNPROTECT(nprotect);
      error("variable %s not specified",names[k]);
    }
    idx[k] -= 1;
  }
  UNPROTECT(nprotect);
  return index;
}

static R_INLINE void setrownames (SEXP x, SEXP names, int n) {
  int nprotect = 0;
  SEXP dimnms, nm;
  PROTECT(nm = AS_CHARACTER(names)); nprotect++;
  PROTECT(dimnms = allocVector(VECSXP,n)); nprotect++;
  SET_ELEMENT(dimnms,0,nm);	// set row names
  SET_DIMNAMES(x,dimnms);
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

static R_INLINE double expit (double x) {
  return 1.0/(1.0 + exp(-x));
}

static R_INLINE double logit (double x) {
  return log(x/(1-x));
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
  while (list != R_NilValue) {
    if (strcmp(STRING_PTR(PRINTNAME(TAG(list))),name)==0) break;
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
