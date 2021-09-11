// -*- C++ -*-

#ifndef _POMP_DEFINES_H_
#define _POMP_DEFINES_H_

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

#include "pomp.h"

# define MATCHROWNAMES(X,N,W) (matchnames(GET_ROWNAMES(GET_DIMNAMES(X)),(N),(W)))
# define MATCHCOLNAMES(X,N,W) (matchnames(GET_COLNAMES(GET_DIMNAMES(X)),(N),(W)))

typedef enum {undef=0,Rfun=1,native=2,regNative=3} pompfunmode;
typedef enum {dflt=0,onestep=1,discrete=2,euler=3,gill=4} rprocmode;

// lookup-table structure, as used internally
typedef struct {
  int length, width;
  int index;
  int order;
  double *x;
  double *y;
} lookup_table_t;

typedef SEXP pomp_fun_handler_t (SEXP pfun, SEXP gnsi, pompfunmode *mode, SEXP S, SEXP P, SEXP O, SEXP C);
typedef SEXP load_stack_incr_t (SEXP pack);
typedef SEXP load_stack_decr_t (SEXP pack);
typedef lookup_table_t make_covariate_table_t (SEXP object, int *ncovar);
typedef void table_lookup_t (lookup_table_t *tab, double x, double *y);
typedef SEXP apply_probe_data_t (SEXP object, SEXP probes);
typedef SEXP apply_probe_sim_t (SEXP object, SEXP nsim, SEXP params, SEXP probes, SEXP datval, SEXP gnsi);
typedef SEXP systematic_resampling_t (SEXP weights);
typedef void set_pomp_userdata_t (SEXP userdata);
typedef void unset_pomp_userdata_t (void);
typedef SEXP get_covariate_names_t (SEXP object);

static R_INLINE SEXP makearray (int rank, const int *dim) {
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

// check if names exist and are nonempty
static R_INLINE int invalid_names (SEXP names) {
  int i, ok;
  ok = !isNull(names);
  for (i = 0; ok && i < LENGTH(names); i++)
    ok = ok && LENGTH(STRING_ELT(names,i)) > 0 &&
      STRING_ELT(names,i) != NA_STRING;
  return !ok;
}

static R_INLINE SEXP matchnames (SEXP provided, SEXP needed, const char *where) {
  int m = LENGTH(provided);
  int n = length(needed);
  SEXP index;
  int *idx, i, j;

  PROTECT(provided = AS_CHARACTER(provided));
  PROTECT(needed = AS_CHARACTER(needed));
  if (invalid_names(provided)) errorcall(R_NilValue,"invalid variable names among the %s.",where);
  PROTECT(index = NEW_INTEGER(n));
  idx = INTEGER(index);

  for (i = 0; i < n; i++) {
    for (j = 0; j < m; j++) {
      if (!strcmp(CHAR(STRING_ELT(provided,j)),CHAR(STRING_ELT(needed,i)))) {
        idx[i] = j;
        break;
      }
    }
    if (j==m) errorcall(R_NilValue,"variable '%s' not found among the %s.",CHAR(STRING_ELT(needed,i)),where);
  }
  UNPROTECT(3);
  return index;
}

static R_INLINE void setrownames (SEXP x, SEXP names, int rank) {
  SEXP dimnms, nm;
  PROTECT(nm = AS_CHARACTER(names));
  PROTECT(dimnms = allocVector(VECSXP,rank));
  SET_ELEMENT(dimnms,0,nm);	// set row names
  SET_DIMNAMES(x,dimnms);
  UNPROTECT(2);
}

// This only works if the dimnames have already been created and set
// e.g., with 'setrownames'
static R_INLINE void setcolnames (SEXP x, SEXP names) {
  SEXP dn;
  PROTECT(dn = GET_DIMNAMES(x));
  SET_ELEMENT(dn,1,names);
  SET_DIMNAMES(x,dn);
  UNPROTECT(1);
}

static R_INLINE void fixdimnames (SEXP x, const char **names, int n) {
  int nprotect = 2;
  int i;
  SEXP dimnames, nm;
  PROTECT(dimnames = GET_DIMNAMES(x));
  if (isNull(dimnames)) {
    PROTECT(dimnames = allocVector(VECSXP,n)); nprotect++;
  }
  PROTECT(nm = allocVector(VECSXP,n));
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
    if (strcmp(CHAR(STRING_ELT(names,i)),str) == 0) {
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

static R_INLINE SEXP paste0 (SEXP a, SEXP b) {
  SEXP x;
  PROTECT(x = LCONS(b,R_NilValue));
  PROTECT(x = LCONS(a,x));
  PROTECT(x = LCONS(install("paste0"),x));
  PROTECT(x = eval(x,R_BaseEnv));
  UNPROTECT(4);
  return x;
}

static R_INLINE SEXP paste (SEXP a, SEXP b, SEXP sep) {
  SEXP x;
  PROTECT(x = LCONS(sep,R_NilValue));
  SET_TAG(x,install("sep"));
  PROTECT(x = LCONS(b,x));
  PROTECT(x = LCONS(a,x));
  PROTECT(x = LCONS(install("paste"),x));
  PROTECT(x = eval(x,R_BaseEnv));
  UNPROTECT(5);
  return x;
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
