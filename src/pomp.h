// -*- C++ -*-
// Header for the C API for pomp.
// Documentation: https://kingaa.github.io/pomp/vignettes/C_API.html

#ifndef _POMP_H_
#define _POMP_H_

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

typedef void periodic_bspline_basis_eval_t (double x, double period, int degree, int nbasis, double *y);
typedef void periodic_bspline_basis_eval_deriv_t (double x, double period, int degree, int nbasis, int deriv, double *y);
typedef const SEXP get_userdata_t (const char *name);
typedef const int *get_userdata_int_t (const char *name);
typedef const double *get_userdata_double_t (const char *name);

static R_INLINE double rgammawn (double sigma, double dt) {
  double sigmasq;
  sigmasq = sigma*sigma;
  return (sigmasq > 0) ? rgamma(dt/sigmasq,sigmasq) : dt;
}

static R_INLINE double dot_product (int dim, const double *basis, const double *coef) {
  int j;
  double trans = 0.0;
  for (j = 0; j < dim; j++)
    trans += coef[j]*basis[j];
  return(trans);
}


static R_INLINE double logit (double p) {
  return log(p/(1.0-p));
}

static R_INLINE double expit (double x) {
  return 1.0/(1.0+exp(-x));
}

static R_INLINE void reulermultinom (int m, double size, const double *rate,
  double dt, double *trans) {
  double p = 0.0;
  int j, k;
  if ((size < 0.0) || (dt < 0.0) || (floor(size+0.5) != size)) {
    for (k = 0; k < m; k++) trans[k] = R_NaN;
    return;
  }
  for (k = 0; k < m; k++) {
    if (rate[k] < 0.0) {
      for (j = 0; j < m; j++) trans[j] = R_NaN;
      return;
    }
    p += rate[k]; // total event rate
  }
  if (p > 0.0) {
    size = rbinom(size,1-exp(-p*dt)); // total number of events
    if (!(R_FINITE(size)))
      warningcall(R_NilValue,"in 'reulermultinom': result of binomial draw is not finite.");
    m -= 1;
    for (k = 0; k < m; k++) {
      if (rate[k] > p) p = rate[k];
      trans[k] = ((size > 0) && (p > 0)) ? rbinom(size,rate[k]/p) : 0;
      if (!(R_FINITE(size)&&R_FINITE(p)&&R_FINITE(rate[k])&&R_FINITE(trans[k])))
        warningcall(R_NilValue,"in 'reulermultinom': result of binomial draw is not finite.");
      size -= trans[k];
      p -= rate[k];
    }
    trans[m] = size;
  } else {
    for (k = 0; k < m; k++) trans[k] = 0.0;
  }
}

static R_INLINE double deulermultinom (int m, double size, const double *rate,
  double dt, double *trans, int give_log) {
  double p = 0.0;
  double n = 0.0;
  double ff = 0.0;
  int k;
  if ((dt < 0.0) || (size < 0.0) || (floor(size+0.5) != size)) {
    warningcall(R_NilValue,"in 'deulermultinom': NaNs produced.");
    return R_NaN;
  }
  for (k = 0; k < m; k++) {
    if (rate[k] < 0.0) {
      warningcall(R_NilValue,"in 'deulermultinom': NaNs produced.");
      return R_NaN;
    }
    if (trans[k] < 0.0) {
      ff = (give_log) ? R_NegInf: 0.0;
      return ff;
    }
    p += rate[k]; // total event rate
    n += trans[k]; // total number of events
  }
  if (n > size) {
    ff = (give_log) ? R_NegInf: 0.0;
    return ff;
  }
  ff = dbinom(n,size,1-exp(-p*dt),1); // total number of events
  m -= 1;
  for (k = 0; k < m; k++) {
    if ((n > 0) && (p > 0)) {
      if (rate[k] > p) p = rate[k];
      ff += dbinom(trans[k],n,rate[k]/p,1);
    }
    n -= trans[k];
    p -= rate[k];
  }
  ff = (give_log) ? ff : exp(ff);
  return ff;
}

static R_INLINE double dmultinom (int m, const double *prob, double *x, int give_log) {
  double p = 0.0;
  double n = 0.0;
  double ff = 0.0;
  int k;

  for (k = 0; k < m; k++) {
    if (prob[k] < 0.0) {
      warningcall(R_NilValue,"in 'dmultinom': NaNs produced.");
      return R_NaN;
    }

    if ((x[k] < 0.0) || (floor(x[k]+0.5) != x[k])) {
      ff = (give_log) ? R_NegInf: 0.0;
      return ff;
    }

    p += prob[k]; // sum of probabilities
    n += x[k]; // total number of events
  }

  for (k = 0; k < m; k++) {
    if ((n > 0) && (p > 0)) {
      if (prob[k] > p) p = prob[k];
      ff += dbinom(x[k],n,prob[k]/p,1);
    } else if (x[k] < 0.0) {
      ff = R_NegInf;
      return ff;
    }

    n -= x[k];
    p -= prob[k];
  }

  ff = (give_log) ? ff : exp(ff);
  return ff;
}

static R_INLINE void to_log_barycentric (double *xt, const double *x, int n) {
  double sum;
  int i;
  for (i = 0, sum = 0.0; i < n; i++) sum += x[i];
  for (i = 0; i < n; i++) xt[i] = log(x[i]/sum);
}

static R_INLINE void from_log_barycentric (double *xt, const double *x, int n) {
  double sum;
  int i;
  for (i = 0, sum = 0.0; i < n; i++) sum += (xt[i] = exp(x[i]));
  for (i = 0; i < n; i++) xt[i] /= sum;
}

static R_INLINE double exp2geom_rate_correction (double R, double dt) {
  return (dt > 0) ? log(1.0+R*dt)/dt : R;
}

static R_INLINE double rbetabinom (double size, double prob, double theta) {
  return rbinom(size,rbeta(prob*theta,(1.0-prob)*theta));
}

static R_INLINE double dbetabinom (double x, double size, double prob,
  double theta, int give_log) {
  double a = theta*prob;
  double b = theta*(1.0-prob);
  double f = lchoose(size,x)-lbeta(a,b)+lbeta(a+x,b+size-x);
  return (give_log) ? f : exp(f);
}

static R_INLINE double rbetanbinom (double mu, double size, double theta) {
  double p = size/(size+mu);
  return rnbinom(size,rbeta(p*theta,(1.0-p)*theta));
}

static R_INLINE double dbetanbinom (double x, double mu, double size,
  double theta, int give_log) {
  double p = size/(size+mu);
  double a = theta*p;
  double b = theta*(1.0-p);
  double f = lchoose(size+x-1,size-1)-lbeta(a,b)+lbeta(a+size,b+x);
  return (give_log) ? f : exp(f);
}

typedef void pomp_rinit(double *x, const double *p, double t,
  const int *stateindex, const int *parindex, const int *covindex,
  const double *covars);

typedef double pomp_ssa_rate_fn(int event, double t, const double *x, const double *p,
  const int *stateindex, const int *parindex, const int *covindex, const double *covars);

typedef void pomp_onestep_sim(double *x, const double *p,
  const int *stateindex, const int *parindex, const int *covindex,
  const double *covars, double t, double dt);

typedef void pomp_onestep_pdf(double *loglik,
  const double *x1, const double *x2, double t1, double t2, const double *p,
  const int *stateindex, const int *parindex, const int *covindex,
  const double *covars);

typedef void pomp_skeleton (double *f, const double *x, const double *p,
  const int *stateindex, const int *parindex, const int *covindex,
  const double *covars, double t);

typedef void pomp_measure_model_simulator (double *y, const double *x, const double *p,
  const int *obsindex, const int *stateindex, const int *parindex, const int *covindex,
  const double *covars, double t);

typedef void pomp_measure_model_expectation (double *f, const double *x, const double *p,
  const int *obsindex, const int *stateindex, const int *parindex, const int *covindex,
  const double *covars, double t);

typedef void pomp_measure_model_covariance (double *f, const double *x, const double *p,
  const int *vmatindex, const int *stateindex, const int *parindex, const int *covindex,
  const double *covars, double t);

typedef void pomp_measure_model_density (double *lik, const double *y, const double *x, const double *p, int give_log,
					 const int *obsindex, const int *stateindex, const int *parindex, const int *covindex,
					 const double *covars, double t);

typedef void pomp_rprior (double *p, const int *parindex);

typedef void pomp_dprior (double *lik, const double *p, int give_log, const int *parindex);

typedef void pomp_transform_fn (double *pt, const double *p, const int *parindex);

#endif
