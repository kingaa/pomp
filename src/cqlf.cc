#include <R.h>
#include <Rmath.h>
#include <float.h>
#include "pomp.h"

static inline int max (int a, int b) {
  return((a > b) ? a : b);
}

extern "C" {

  void composite_quasilikelihood (double *Q, 
				  int *dim, int *ldim,
				  double *x, double *y, 
				  int *maxlag
				  )
  {
    int T, N, B;

    if (*ldim > 2) {
      T = dim[0];
      N = dim[1];
      B = dim[2];
    } else {
      T = dim[0];
      N = 1;
      B = dim[1];
    }

    view<double> X(x,T,N), Y(y,T);
    double logss[T][N];

    for (int j = 0; j < N; j++) {
      for (int t = 0; t < T; t++) {
	// compute the mean of the stochastic realizations
	double mu, sigma, tmp, sum = 0.0;
	for (int b = 0; b < B; b++) sum += X(t,j,b);
	mu = sum / ((double) B);
	// compute the SD of the stochastic realizations
	// center the stochastic realizations as we go
	sum = 0.0;
	for (int b = 0; b < B; b++) {
	  X(t,j,b) -= mu;
	  tmp = X(t,j,b);
	  sum += tmp * tmp;
	}
	sigma = sqrt(sum / ((double) B)) + DBL_EPSILON;
	logss[t][j] = 2*log(sigma);

	// standardize the stochastic realizations
	for (int b = 0; b < B; b++) X(t,j,b) /= sigma;
      
	// standardize the data
	Y(t,j) -= mu;
	Y(t,j) /= sigma;
      }
    }

    // compute the quasilikelihood
    *Q = 0;
    for (int j = 0; j < N; j++) {
      for (int t = 0; t < T; t++) {
	for (int s = max(t-maxlag[0],0); s < t; s++) {
	  double tmp, rho = 0.0;
	  for (int b = 0; b < B; b++) rho += X(s,j,b)*X(t,j,b);
	  rho /= ((double) B);
	  tmp = fabs(1.0 - rho*rho) + DBL_EPSILON;
	  *Q += logss[s][j] + logss[t][j] + log(tmp) 
	    + (Y(s,j)*Y(s,j) + Y(t,j)*Y(t,j) - 2*rho*Y(s,j)*Y(t,j))/tmp;
	}
      } 
      for (int k = 0; k < N; k++) {
	if (k != j) {
	  for (int t = 0; t < T; t++) {
	    for (int s = max(t-maxlag[1],0); s <= t; s++) {
	      if ((s != t) || (k < j)) {
		double tmp, rho = 0.0;
		for (int b = 0; b < B; b++) rho += X(s,k,b)*X(t,j,b);
		rho /= ((double) B);
		tmp = fabs(1.0 - rho*rho) + DBL_EPSILON;
		*Q += logss[s][k] + logss[t][j] + log(tmp) 
		  + (Y(s,k)*Y(s,k) + Y(t,j)*Y(t,j) - 2*rho*Y(s,k)*Y(t,j))/tmp;
	      }
	    }
	  }
	}
      }
    }
  }

}
