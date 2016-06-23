// -*- mode: C++ -*-

#include "pomp_internal.h"
#include "pomp_mat.h"
#include <stdio.h>

static void order_reg_solve (double *beta, double *x, double *mm, double *tau, int *pivot, int n, int np, int diff);
static void order_reg_model_matrix (double *z, double *X, double *tau, int *pivot, int n, int np, int diff);

SEXP probe_marginal_setup (SEXP ref, SEXP order, SEXP diff) {
  int nprotect = 0;
  SEXP z, mm, tau, pivot, retval, retvalnames;
  int n, nx, np, df, dim[2];

  np = *(INTEGER(AS_INTEGER(order))); // order of polynomial regression
  df = *(INTEGER(AS_INTEGER(diff)));  // order of differencing
  n = LENGTH(ref);		      // n = number of observations
  nx = n - df;			      // nx = rows in model matrix
  dim[0] = nx; dim[1] = np;	      // dimensions of model matrix

  if (nx < 1) errorcall(R_NilValue,"must have diff < number of observations");

  PROTECT(z = duplicate(AS_NUMERIC(ref))); nprotect++;
  PROTECT(mm = makearray(2,dim)); nprotect++;
  PROTECT(tau = NEW_NUMERIC(np)); nprotect++;
  PROTECT(pivot = NEW_INTEGER(np)); nprotect++;

  PROTECT(retval = NEW_LIST(3)); nprotect++;
  PROTECT(retvalnames = NEW_CHARACTER(3)); nprotect++;
  SET_STRING_ELT(retvalnames,0,mkChar("mm"));
  SET_STRING_ELT(retvalnames,1,mkChar("tau"));
  SET_STRING_ELT(retvalnames,2,mkChar("pivot"));
  SET_ELEMENT(retval,0,mm);
  SET_ELEMENT(retval,1,tau);
  SET_ELEMENT(retval,2,pivot);
  SET_NAMES(retval,retvalnames);

  order_reg_model_matrix(REAL(z),REAL(mm),REAL(tau),INTEGER(pivot),n,np,df);
  
  UNPROTECT(nprotect);
  return(retval);
}

SEXP probe_marginal_solve (SEXP x, SEXP setup, SEXP diff) {
  int nprotect = 0;
  SEXP X, mm, tau, pivot, beta, beta_names;
  int n, nx, np, df;
  int i;
  char tmp[BUFSIZ];

  df = *(INTEGER(AS_INTEGER(diff)));  // order of differencing
  n = LENGTH(x);		      // n = number of observations

  // unpack the setup information
  PROTECT(mm = VECTOR_ELT(setup,0)); nprotect++; //  QR-decomposed model matrix
  PROTECT(tau = VECTOR_ELT(setup,1)); nprotect++; // diagonals
  PROTECT(pivot = VECTOR_ELT(setup,2)); nprotect++; // pivots

  nx = INTEGER(GET_DIM(mm))[0];	// nx = number of rows in model matrix
  np = INTEGER(GET_DIM(mm))[1];	// np = order of polynomial
  
  if (n-df != nx) errorcall(R_NilValue,"length of 'ref' must equal length of data");
  PROTECT(X = duplicate(AS_NUMERIC(x))); nprotect++; 
   
  PROTECT(beta = NEW_NUMERIC(np)); nprotect++;
  PROTECT(beta_names = NEW_STRING(np)); nprotect++;
  for (i = 0; i < np; i++) {
    snprintf(tmp,BUFSIZ,"marg.%d",i+1);
    SET_STRING_ELT(beta_names,i,mkChar(tmp));
  }
  SET_NAMES(beta,beta_names);

  order_reg_solve(REAL(beta),REAL(X),REAL(mm),REAL(tau),INTEGER(pivot),n,np,df);

  UNPROTECT(nprotect);
  return(beta);
}

// thanks to Simon N. Wood for the original version of the following code
static void order_reg_model_matrix (double *z, double *X, double *tau, int *pivot, int n, int np, int diff) {
  //   z is an n vector, containing no NAs
  //   X is an (n-diff) by np double matrix
  //   pivot is an (n-diff) integer vector
  //   tau is an (n-diff) double vector
  //   This routine first differences z 'diff' times.
  //   z is centred and sorted
  //   then the model matrix is set up and QR decomposed
  int nx;
  double *X1, *X2, xx;
  int i, j;

  for (i = 0, nx = n; i < diff; i++) { // differencing loop
    nx--;
    for (j = 0; j < nx; j++) z[j] = z[j+1] - z[j];
  }
  // nx = number of rows in model matrix
  // z is now an nx-vector

  // center z
  for (j = 0, xx = 0.0; j < nx; j++) xx += z[j]; xx /= nx; // xx = mean(z)
  for (j = 0; j < nx; j++) z[j] -= xx;

  // now sort
  R_qsort(z,1,nx);

  // Now create the model matrix, using contents of z 
  if (np < 1) np = 1;
  for (j = 0; j < nx; j++) X[j] = z[j]; // first column
  for (i = 1, X1 = X, X2 = X+nx; i < np; i++, X1 += nx, X2 += nx) 
    for (j = 0; j < nx; j++) X2[j] = X1[j]*z[j];

  // QR decompose the model matrix 
  pomp_qr(X,nx,np,pivot,tau);
  
}

// thanks to Simon N. Wood for the original version of the following code
static void order_reg_solve (double *beta, double *x, double *mm, double *tau, int *pivot, int n, int np, int diff) {
  int nx;
  double xx;
  int i, j;

  for (i = 0, nx = n; i < diff; i++) { // differencing loop
    nx--;
    for (j = 0; j < nx; j++) x[j] = x[j+1] - x[j];
  }
  // nx = number of rows in model matrix
  // x is now an nx-vector

  // center x
  for (j = 0, xx = 0.0; j < nx; j++) xx += x[j]; xx /= nx; // xx = mean(x)
  for (j = 0; j < nx; j++) x[j] -= xx;

  // now sort
  R_qsort(x,1,nx);

  // solve R b = Q'x for b
  pomp_qrqy(x,mm,nx,tau,nx,1,np,"left","transpose"); // y <- Q'y
  pomp_backsolve(mm,nx,np,x,1,"Upper","No transpose","Non-unit");   // y <- R^{-1} y

  // unpivot and store the coefficients in beta
  for (i = 0; i < np; i++) beta[pivot[i]] = x[i];

}
