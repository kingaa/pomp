// -*- C++ -*-

#include <Rdefines.h>
#include <R_ext/Constants.h>

#include "pomp_internal.h"

SEXP iterate_map (SEXP object, SEXP times, SEXP t0, SEXP x0, SEXP params)
{
  int nprotect = 0;
  SEXP fn, fcall, rho, ans, nm;
  SEXP X;
  SEXP Snames, Pnames, Cnames;
  SEXP zeronames;
  SEXP cvec, tvec, xvec, pvec;
  int nvars, npars, nreps, nrepx, nrepp, ntimes, ncovars, nzeros, nsteps;
  int mode = -1;
  int *dim;
  int *zidx = 0, *sidx = 0, *pidx = 0, *cidx = 0;
  pomp_skeleton *ff = NULL;
  struct lookup_table covariate_table;

  PROTECT(x0 = as_matrix(x0)); nprotect++;
  dim = INTEGER(GET_DIM(x0));
  nvars = dim[0]; nrepx = dim[1];

  PROTECT(params = as_matrix(params)); nprotect++;
  dim = INTEGER(GET_DIM(params));
  npars = dim[0]; nrepp = dim[1];

  // 2nd dimension of 'x0' and 'params' need not agree
  nreps = (nrepp > nrepx) ? nrepp : nrepx;
  if ((nreps % nrepp != 0) || (nreps % nrepx != 0))
    error("trajectory error: 2nd dimensions of 'x0' and 'params' are incompatible");

  PROTECT(times = AS_NUMERIC(times)); nprotect++;
  ntimes = LENGTH(times);

  PROTECT(Snames = GET_ROWNAMES(GET_DIMNAMES(x0))); nprotect++;
  PROTECT(Pnames = GET_ROWNAMES(GET_DIMNAMES(params))); nprotect++;
  PROTECT(Cnames = GET_COLNAMES(GET_DIMNAMES(GET_SLOT(object,install("covar"))))); nprotect++;
    
  // set up the covariate table
  covariate_table = make_covariate_table(object,&ncovars);

  // vector for interpolated covariates
  PROTECT(cvec = NEW_NUMERIC(ncovars)); nprotect++;
  SET_NAMES(cvec,Cnames);

  // extract user-defined function
  PROTECT(fn = pomp_fun_handler(GET_SLOT(object,install("skeleton")),&mode)); nprotect++;

  // extract 'userdata' as pairlist
  PROTECT(fcall = VectorToPairList(GET_SLOT(object,install("userdata")))); nprotect++;

  // get names and indices of accumulator variables
  PROTECT(zeronames = GET_SLOT(object,install("zeronames"))); nprotect++;
  nzeros = LENGTH(zeronames);
  if (nzeros > 0) {
    zidx = INTEGER(PROTECT(matchnames(Snames,zeronames))); nprotect++;
  }

  // set up the computations
  switch (mode) {

  case 0:                       // R function
    
    PROTECT(tvec = NEW_NUMERIC(1)); nprotect++;
    PROTECT(xvec = NEW_NUMERIC(nvars)); nprotect++;
    PROTECT(pvec = NEW_NUMERIC(npars)); nprotect++;
    SET_NAMES(xvec,Snames);
    SET_NAMES(pvec,Pnames);

    // set up the function call
    PROTECT(fcall = LCONS(cvec,fcall)); nprotect++;
    SET_TAG(fcall,install("covars"));
    PROTECT(fcall = LCONS(pvec,fcall)); nprotect++;
    SET_TAG(fcall,install("params"));
    PROTECT(fcall = LCONS(tvec,fcall)); nprotect++;
    SET_TAG(fcall,install("t"));
    PROTECT(fcall = LCONS(xvec,fcall)); nprotect++;
    SET_TAG(fcall,install("x"));
    PROTECT(fcall = LCONS(fn,fcall)); nprotect++;

    // get function's environment
    PROTECT(rho = (CLOENV(fn))); nprotect++;

    break;

  case 1:                       // native skeleton
    
    // construct state, parameter, covariate, observable indices
    sidx = INTEGER(PROTECT(name_index(Snames,object,"statenames"))); nprotect++;
    pidx = INTEGER(PROTECT(name_index(Pnames,object,"paramnames"))); nprotect++;
    cidx = INTEGER(PROTECT(name_index(Cnames,object,"covarnames"))); nprotect++;
    
    ff = (pomp_skeleton *) R_ExternalPtrAddr(fn);
    
    break;

  default:
    error("unrecognized 'mode' in 'iterate_map'");
    break;
  }

  // create array to store results
  {
    int dim[3] = {nvars, nreps, ntimes};
    PROTECT(X = makearray(3,dim)); nprotect++;
    setrownames(X,Snames,3);
  }

  switch (mode) {

  case 0:			// R function

    {
      int first = 1;
      int use_names;
      double *time = REAL(times);
      double *x0s = REAL(x0);
      double *Xt = REAL(X);
      double *ps = REAL(params);
      double *cp = REAL(cvec);
      double *tp = REAL(tvec);
      double *xp = REAL(xvec);
      double *pp = REAL(pvec);
      double *fs, *xm;
      double deltat = *(REAL(GET_SLOT(object,install("skelmap.delta.t"))));
      double t = *(REAL(AS_NUMERIC(t0)));
      int *posn = 0;
      int h, i, j, k;

      for (j = 0; j < nreps; j++)
	for (i = 0; i < nvars; i++) 
	  Xt[i+nvars*j] = x0s[i+nvars*(j%nrepx)];

      for (k = 0; k < ntimes; k++, time++, Xt += nvars*nreps) {
	
	R_CheckUserInterrupt();
	    
	nsteps = num_map_steps(t,*time,deltat); 

	for (i = 0; i < nzeros; i++)
	  for (j = 0, xm = &Xt[zidx[i]]; j < nreps; j++, xm += nvars)
	    *xm = 0.0;
	
	for (h = 0; h < nsteps; h++) {
	  
	  table_lookup(&covariate_table,t,cp,0);
	  
	  for (j = 0, xm = Xt; j < nreps; j++, xm += nvars) {
	    
	    *tp = t;
	    memcpy(xp,xm,nvars*sizeof(double));
	    memcpy(pp,&ps[npars*(j%nrepp)],npars*sizeof(double));
	    
	    if (first) {
	      
	      PROTECT(ans = eval(fcall,rho)); nprotect++;
	      if (LENGTH(ans)!=nvars)
		error("user 'skeleton' returns a vector of %d state variables but %d are expected",LENGTH(ans),nvars);

	      // get name information to fix possible alignment problems
	      PROTECT(nm = GET_NAMES(ans)); nprotect++;
	      use_names = !isNull(nm);
	      if (use_names) {
		posn = INTEGER(PROTECT(matchnames(Snames,nm))); nprotect++;
	      }

	      fs = REAL(AS_NUMERIC(ans));
          
	      first = 0;
	      
	    } else {
	      
	      fs = REAL(AS_NUMERIC(eval(fcall,rho)));
	      
	    }

	    if (use_names) 
	      for (i = 0; i < nvars; i++) xm[posn[i]] = fs[i];
	    else
	      for (i = 0; i < nvars; i++) xm[i] = fs[i];
 
	  }

	  t += deltat;

	}

	if (k+1 < ntimes)
	  memcpy(Xt+nvars*nreps,Xt,nvars*nreps*sizeof(double));
	
      }
    }
    break;

  case 1: 			// native code

    {
      double *time = REAL(times);
      double *Xt = REAL(X);
      double *ps = REAL(params);
      double *x0s = REAL(x0);
      double *cp = REAL(cvec);
      double deltat = *(REAL(GET_SLOT(object,install("skelmap.delta.t"))));
      double t = *(REAL(AS_NUMERIC(t0)));
      double *xs = (double *) R_alloc(nvars*nreps,sizeof(double)); // array to temporarily hold states
      double *fp, *xp;
      int h, i, j, k;

      // initialize the state matrix
      for (j = 0; j < nreps; j++)
	for (i = 0; i < nvars; i++) 
	  xs[i+nvars*j] = x0s[i+nvars*(j%nrepx)];

      for (k = 0; k < ntimes; k++, time++, Xt += nvars*nreps) {
	
	R_CheckUserInterrupt();

	for (i = 0; i < nzeros; i++)
	  for (j = 0, xp = &Xt[zidx[i]]; j < nreps; j++, xp += nvars)
	    *xp = 0.0;
	
	nsteps = num_map_steps(t,*time,deltat);

	for (h = 0; h < nsteps; h++) {
	  
	  table_lookup(&covariate_table,t,cp,0);
	  
	  for (j = 0, fp = Xt, xp = xs; j < nreps; j++, fp += nvars, xp += nvars) {
	    (*ff)(fp,xp,&ps[npars*(j%nrepp)],sidx,pidx,cidx,ncovars,cp,t);
	    memcpy(xp,fp,nvars*sizeof(double));
	  }

	  t += deltat;

	}

	memcpy(Xt,xs,nvars*nreps*sizeof(double));

      }
    }
    break;

  default:
    error("unrecognized 'mode' in 'iterate_map'");
    break;
  }

  UNPROTECT(nprotect);
  return X;
}

static struct {
  SEXP object;
  SEXP params;
  SEXP skelfun;
  SEXP xnames;
  int xdim[3];
} _pomp_vf_eval_common;


#define COMMON(X)    (_pomp_vf_eval_common.X)

void pomp_desolve_setup (SEXP object, SEXP params, SEXP fun, SEXP statenames, SEXP nvar, SEXP nrep) {
  COMMON(object) = object;
  COMMON(params) = params;
  COMMON(skelfun) = fun;
  COMMON(xnames) = statenames;
  COMMON(xdim)[0] = *INTEGER(AS_INTEGER(nvar));
  COMMON(xdim)[1] = *INTEGER(AS_INTEGER(nrep));
  COMMON(xdim)[2] = 1;
}

void pomp_desolve_takedown (void) {
  COMMON(object) = R_NilValue;
  COMMON(params) = R_NilValue;
  COMMON(skelfun) = R_NilValue;
  COMMON(xnames) = R_NilValue;
  COMMON(xdim)[0] = 0;
  COMMON(xdim)[1] = 0;
  COMMON(xdim)[2] = 0;
}

void pomp_vf_eval (int *neq, double *t, double *y, double *ydot, double *yout, int *ip) 
{
  SEXP T, X, dXdt;
  
  PROTECT(T = NEW_NUMERIC(1));
  PROTECT(X = makearray(3,COMMON(xdim)));
  setrownames(X,COMMON(xnames),3);
  REAL(T)[0] = *t;
  memcpy(REAL(X),y,(*neq)*sizeof(double));
  PROTECT(dXdt = do_skeleton(COMMON(object),X,T,COMMON(params),COMMON(skelfun)));
  memcpy(ydot,REAL(dXdt),(*neq)*sizeof(double));
  
  UNPROTECT(3);
}

#undef COMMON

// copy t(x[-1,-1]) -> y[,rep,]
SEXP traj_transp_and_copy (SEXP y, SEXP x, SEXP rep) {
  int nprotect = 0;
  SEXP ans = R_NilValue;
  int nvars, nreps, ntimes;
  int i, j, k, m, n;
  double *xp, *yp;

  j = INTEGER(rep)[0]-1;
  nvars = INTEGER(GET_DIM(y))[0];
  nreps = INTEGER(GET_DIM(y))[1];
  ntimes = INTEGER(GET_DIM(y))[2];
  n = INTEGER(GET_DIM(x))[0];
  m = nvars*nreps;

  for (i = 0, xp = REAL(x)+n+1; i < nvars; i++, xp += n) {
    for (k = 0, yp = REAL(y)+i+nvars*j; k < ntimes; k++, yp += m) {
      *yp = xp[k];
    }
  }

  UNPROTECT(nprotect);
  return ans;
}

