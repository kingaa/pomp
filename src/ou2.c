// -*- C++ -*-

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>

struct lookup_table {
  int length, width;
  int index;
  double *x;
  double **y;
};

// prototypes
SEXP ou2_simulator (SEXP xstart, SEXP times, SEXP params);
SEXP ou2_density (SEXP x, SEXP times, SEXP params, SEXP give_log);
SEXP bivariate_normal_rmeasure (SEXP x, SEXP times, SEXP params);
SEXP bivariate_normal_dmeasure (SEXP y, SEXP x, SEXP times, SEXP params, SEXP give_log);

static double expit (double x);
static double logit (double x);
static void normal_rmeasure (int *n, double *X, double *par, int *index, double *obs);
static void normal_dmeasure (int *n, double *X, double *par, int *index, double *Y, double *f, int *give_log);
static void ou2_adv (double *x, double *xstart, double *par, double *times, int *n, int *parindex);
static void ou2_pdf (double *d, double *X, double *par, double *times, int *n, int *parindex, int *give_log);
static void sim_ou2 (double *x,
		     double alpha1, double alpha2, double alpha3, double alpha4, 
		     double sigma1, double sigma2, double sigma3);
static double dens_ou2 (double *x1, double *x2,
			double alpha1, double alpha2, double alpha3, double alpha4, 
			double sigma1, double sigma2, double sigma3, int give_log);
static SEXP makearray (int rank, int *dim);
static void setrownames (SEXP x, int n, char **names);
static SEXP matchrownames (SEXP x, int n, char **names);

// this is the rprocess function
// it is basically a wrapper around a call to 'ou2_adv', which could be called from R directly
SEXP ou2_simulator (SEXP xstart, SEXP times, SEXP params) {
  int nprotect = 0;
  int *dim, xdim[3], ndim[4];
  int nvar, npar, nrep, ntimes;
  char *paramnames[] = {"alpha.1","alpha.2","alpha.3","alpha.4","sigma.1","sigma.2","sigma.3"};
  int nparams = sizeof(paramnames)/sizeof(paramnames[0]);
  SEXP X, pindex;
  dim = INTEGER(GET_DIM(xstart)); nvar = dim[0]; nrep = dim[1];
  dim = INTEGER(GET_DIM(params)); npar = dim[0];
  ntimes = length(times);
  xdim[0] = nvar; xdim[1] = nrep; xdim[2] = ntimes;
  PROTECT(X = makearray(3,xdim)); nprotect++;
  PROTECT(pindex = matchrownames(params,nparams,paramnames)); nprotect++;
  ndim[0] = nvar; ndim[1] = npar; ndim[2] = nrep; ndim[3] = ntimes;
  ou2_adv(REAL(X),REAL(xstart),REAL(params),REAL(times),ndim,INTEGER(pindex));
  UNPROTECT(nprotect);
  return X;
}

// this is the dprocess function
// it is a wrapper around a call to 'ou2_pdf', which could be called from R directly
SEXP ou2_density (SEXP x, SEXP times, SEXP params, SEXP give_log) {
  int nprotect = 0;
  int *dim, xdim[2], ndim[4];
  int nvar, npar, nrep, ntimes;
  char *paramnames[] = {"alpha.1","alpha.2","alpha.3","alpha.4","sigma.1","sigma.2","sigma.3"};
  int nparams = sizeof(paramnames)/sizeof(paramnames[0]);
  SEXP D, index;
  dim = INTEGER(GET_DIM(x)); nvar = dim[0]; nrep = dim[1];
  dim = INTEGER(GET_DIM(params)); npar = dim[0]; 
  ntimes = length(times);
  xdim[0] = nrep; xdim[1] = ntimes-1;
  PROTECT(D = makearray(2,xdim)); nprotect++;
  PROTECT(index = matchrownames(params,nparams,paramnames)); nprotect++;
  ndim[0] = nvar; ndim[1] = npar; ndim[2] = nrep; ndim[3] = ntimes;
  ou2_pdf(REAL(D),REAL(x),REAL(params),REAL(times),ndim,INTEGER(index),LOGICAL(give_log));
  UNPROTECT(nprotect);
  return D;
}

// this is the rmeasure function
// it is a wrapper around a call to 'normal_rmeasure', which could be called from R directly
SEXP bivariate_normal_rmeasure (SEXP x, SEXP times, SEXP params) {
  int nprotect = 0;
  int *dim, xdim[3], ndim[5];
  int nvar, npar, nrep, ntimes;
  char *paramnames[] = {"tau"}, *obsnames[] = {"y1","y2"};
  int nparams = sizeof(paramnames)/sizeof(paramnames[0]);
  int nobs = sizeof(obsnames)/sizeof(obsnames[0]);
  SEXP obs, index;
  dim = INTEGER(GET_DIM(x)); nvar = dim[0]; nrep = dim[1];
  dim = INTEGER(GET_DIM(params)); npar = dim[0];
  ntimes = length(times);
  xdim[0] = nobs; xdim[1] = nrep; xdim[2] = ntimes;
  PROTECT(obs = makearray(3,xdim)); nprotect++;
  setrownames(obs,nobs,obsnames);
  PROTECT(index = matchrownames(params,nparams,paramnames)); nprotect++;
  ndim[0] = nvar; ndim[1] = npar; ndim[2] = nrep; ndim[3] = ntimes; ndim[4] = nobs;
  normal_rmeasure(ndim,REAL(x),REAL(params),INTEGER(index),REAL(obs));
  UNPROTECT(nprotect);
  return obs;
}

// this is the dmeasure function
// it is a wrapper around a call to 'normal_dmeasure', which could be called from R directly
SEXP bivariate_normal_dmeasure (SEXP y, SEXP x, SEXP times, SEXP params, SEXP give_log) {
  int nprotect = 0;
  int *dim, xdim[2], ndim[5];
  int nobs, nvar, npar, nrep, ntimes;
  char *paramnames[] = {"tau"};
  int nparams = sizeof(paramnames)/sizeof(paramnames[0]);
  SEXP d, index;
  dim = INTEGER(GET_DIM(y)); nobs = dim[0];
  dim = INTEGER(GET_DIM(x)); nvar = dim[0]; nrep = dim[1];
  dim = INTEGER(GET_DIM(params)); npar = dim[0];
  ntimes = length(times);
  xdim[0] = nrep; xdim[1] = ntimes;
  PROTECT(d = makearray(2,xdim)); nprotect++;
  PROTECT(index = matchrownames(params,nparams,paramnames)); nprotect++;
  ndim[0] = nvar; ndim[1] = npar; ndim[2] = nrep; ndim[3] = ntimes; ndim[4] = nobs;
  normal_dmeasure(ndim,REAL(x),REAL(params),INTEGER(index),REAL(y),REAL(d),LOGICAL(give_log));
  UNPROTECT(nprotect);
  return d;
}

#define ALPHA1     (pp[parindex[0]])
#define ALPHA2     (pp[parindex[1]])
#define ALPHA3     (pp[parindex[2]])
#define ALPHA4     (pp[parindex[3]])
#define SIGMA1     (pp[parindex[4]])
#define SIGMA2     (pp[parindex[5]])
#define SIGMA3     (pp[parindex[6]])

// advance the matrix of particles from times[0] to the other times given
// it is assumed that the times are consecutive (FIX THIS!)
static void ou2_adv (double *x, double *xstart, double *par, double *times, int *n, int *parindex)
{
  int nvar = n[0], npar = n[1], nrep = n[2], ntimes = n[3];
  double *xp, *pp;
  int i, j, k;
  GetRNGstate();       // initialize R's pseudorandom number generator
  for (j = 0; j < nrep; j++) {
    xp = &x[nvar*j];		// get address of j-th state vector
    for (i = 0; i < nvar; i++) xp[i] = xstart[i+nvar*j]; // copy xstart into the first slice of x
  }
  for (k = 1; k < ntimes; k++) {
    R_CheckUserInterrupt();
    for (j = 0; j < nrep; j++) {
      xp = &x[nvar*(j+nrep*k)];
      pp = &par[npar*j];
      for (i = 0; i < nvar; i++) xp[i] = x[i+nvar*(j+nrep*(k-1))];
      sim_ou2(xp,ALPHA1,ALPHA2,ALPHA3,ALPHA4,SIGMA1,SIGMA2,SIGMA3); // advance particle
    }
  }
  PutRNGstate();	  // finished with R's random number generator
}

// pdf of a single 2D OU transition
static void ou2_pdf (double *d, double *X, double *par, double *times, int *n, int *parindex, int *give_log)
{
  int nvar = n[0], npar = n[1], nrep = n[2], ntimes = n[3];
  double *x1, *x2, *pp;
  int j, k;
  for (k = 0; k < nrep; k++) {
    pp = &par[npar*k];	       // get address of k-th parameter vector
    x1 = &X[nvar*k];	     // get address of (0,0)-th state vector
    for (j = 1; j < ntimes; j++) {
      R_CheckUserInterrupt();
      x2 = &X[nvar*(k+nrep*j)]; // get address of (k,j)-th state vector
      d[k+nrep*(j-1)] = dens_ou2(x1,x2,ALPHA1,ALPHA2,ALPHA3,ALPHA4,SIGMA1,SIGMA2,SIGMA3,*give_log);
      x1 = x2;
    }
  }
}

#undef ALPHA1
#undef ALPHA2
#undef ALPHA3
#undef ALPHA4
#undef SIGMA1
#undef SIGMA2
#undef SIGMA3

#define TAU   (p[index[0]])

// bivariate normal measurement error density
static void normal_dmeasure (int *n, double *X, double *par, int *index, double *Y, double *f, int *give_log) {
  int nvar = n[0], npar = n[1], nrep = n[2], ntimes = n[3], nobs = n[4];
  double *x, *p, *y, v, val;
  double tol = 1.0e-18;	// tol should be less than the tol in the particle filter!
  int j, k;
  for (j = 0; j < ntimes; j++) {
    R_CheckUserInterrupt();
    y = &Y[nobs*j];
    for (k = 0; k < nrep; k++) {
      x = &X[nvar*(k+j*nrep)];
      p = &par[npar*k];
      v = fabs(TAU);
      if (R_FINITE(v)) {
	val = 0.0;
	if (!ISNA(y[0])) val = dnorm(y[0],x[0],v+tol,1);
	if (!ISNA(y[1])) val += dnorm(y[1],x[1],v+tol,1);
// when give_log=TRUE, use the 1st-order Taylor expansion of log(p+tol) = log(p)+log(1+tol/p) ~= log(p)+tol/p
	f[k+j*nrep] = (*give_log) ? val+tol*exp(-val) : exp(val)+tol; 
      } else {
	f[k+j*nrep] = (*give_log) ? log(tol) : tol;
      }
    }
  }
}

// bivariate normal measurement error simulator
static void normal_rmeasure (int *n, double *X, double *par, int *index, double *obs) {
  int nvar = n[0], npar = n[1], nrep = n[2], ntimes = n[3], nobs = n[4];
  double *x, *p, v;
  double tol = 1.0e-18;	// tol should be less than the tol in the particle filter!
  int j, k;
  GetRNGstate();
  for (j = 0; j < ntimes; j++) {
    R_CheckUserInterrupt();
    for (k = 0; k < nrep; k++) {
      x = &X[nvar*(k+j*nrep)];
      p = &par[npar*k];
      v = fabs(TAU);
      if (R_FINITE(v)) {
	obs[nobs*(k+j*nrep)] = rnorm(x[0],v+tol);
	obs[1+nobs*(k+j*nrep)] = rnorm(x[1],v+tol);
      } else {
	obs[nobs*(k+j*nrep)] = R_NaReal;
	obs[1+nobs*(k+j*nrep)] = R_NaReal;
      }
    }
  }
  PutRNGstate();
}

#undef TAU

// simple 2D Ornstein-Uhlenbeck process simulation
static void sim_ou2 (double *x,
		     double alpha1, double alpha2, double alpha3, double alpha4, 
		     double sigma1, double sigma2, double sigma3)
{
  double eps[2], xnew[2];

  if (!(R_FINITE(x[0]))) return;
  if (!(R_FINITE(x[1]))) return;
  if (!(R_FINITE(alpha1))) return;
  if (!(R_FINITE(alpha2))) return;
  if (!(R_FINITE(alpha3))) return;
  if (!(R_FINITE(alpha4))) return;
  if (!(R_FINITE(sigma1))) return;
  if (!(R_FINITE(sigma2))) return;
  if (!(R_FINITE(sigma3))) return;

  eps[0] = rnorm(0,1);
  eps[1] = rnorm(0,1);

  xnew[0] = alpha1*x[0]+alpha3*x[1]+sigma1*eps[0];
  xnew[1] = alpha2*x[0]+alpha4*x[1]+sigma2*eps[0]+sigma3*eps[1];

  x[0] = xnew[0];
  x[1] = xnew[1];
}

// simple 2D Ornstein-Uhlenbeck process density
static double dens_ou2 (double *x1, double *x2,
			double alpha1, double alpha2, double alpha3, double alpha4, 
			double sigma1, double sigma2, double sigma3, int give_log)
{
  double eps[2], val;

  if (!(R_FINITE(x1[0]))) return R_NaReal;
  if (!(R_FINITE(x1[1]))) return R_NaReal;
  if (!(R_FINITE(x2[0]))) return R_NaReal;
  if (!(R_FINITE(x2[1]))) return R_NaReal;
  if (!(R_FINITE(alpha1))) return R_NaReal;
  if (!(R_FINITE(alpha2))) return R_NaReal;
  if (!(R_FINITE(alpha3))) return R_NaReal;
  if (!(R_FINITE(alpha4))) return R_NaReal;
  if (!(R_FINITE(sigma1))) return R_NaReal;
  if (!(R_FINITE(sigma2))) return R_NaReal;
  if (!(R_FINITE(sigma3))) return R_NaReal;

  // compute residuals
  eps[0] = x2[0]-alpha1*x1[0]-alpha3*x1[1];
  eps[1] = x2[1]-alpha2*x1[0]-alpha4*x1[1];

  // backsolve
  eps[0] /= sigma1;
  eps[1] -= sigma2*eps[0];
  eps[1] /= sigma3;

  val = dnorm(eps[0],0.0,1.0,1)+dnorm(eps[1],0.0,1.0,1);
  return ((give_log) ? val : exp(val));
}

static double expit (double x) {
  return 1.0/(1.0 + exp(-x));
}

static double logit (double x) {
  return log(x/(1-x));
}

static SEXP makearray (int rank, int *dim) {
  int nprotect = 0;
  int *dimp, k;
  SEXP dimx, x;
  PROTECT(dimx = NEW_INTEGER(rank)); nprotect++;
  dimp = INTEGER(dimx); 
  for (k = 0; k < rank; k++) dimp[k] = dim[k];
  PROTECT(x = allocArray(REALSXP,dimx)); nprotect++;
  UNPROTECT(nprotect);
  return x;
}

static void setrownames (SEXP x, int n, char **names) {
  int nprotect = 0;
  int k;
  SEXP nm, dimnms;
  PROTECT(nm = NEW_CHARACTER(n)); nprotect++;
  for (k = 0; k < n; k++) 
    SET_ELEMENT(nm,k,mkChar(names[k]));
  PROTECT(dimnms = allocVector(VECSXP,3)); nprotect++;
  SET_ELEMENT(dimnms,0,nm);	// set row names
  SET_DIMNAMES(x,dimnms);
  UNPROTECT(nprotect);
}

static SEXP matchrownames (SEXP x, int n, char **names) {
  int nprotect = 0;
  int *idx, k;
  SEXP index, nm;
  PROTECT(nm = NEW_CHARACTER(n)); nprotect++;
  for (k = 0; k < n; k++) {
    SET_ELEMENT(nm,k,mkChar(names[k]));
  }
  PROTECT(index = match(GET_ROWNAMES(GET_DIMNAMES(x)),nm,0)); nprotect++;
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

