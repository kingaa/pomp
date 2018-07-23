// dear emacs, please treat this as -*- C++ -*-

#include "pomp_internal.h"
#include <R_ext/Constants.h>

// method: 1 = one-step, 2 = fixed step, 3 = Euler

SEXP euler_model_simulator (SEXP func,
                            SEXP xstart, SEXP times, SEXP params,
                            double deltat, int method, SEXP zeronames,
                            SEXP tcovar, SEXP covar, SEXP args, SEXP gnsi)
{
  int nprotect = 0;
  pompfunmode mode = undef;
  int nvars, npars, nreps, ntimes, nzeros, ncovars, covlen;
  int nstep = 0;
  double dt;
  SEXP X;
  SEXP ans, nm, fn, fcall = R_NilValue, rho = R_NilValue;
  SEXP Snames, Pnames, Cnames;
  SEXP cvec, tvec = R_NilValue;
  SEXP xvec = R_NilValue, pvec = R_NilValue, dtvec = R_NilValue;
  int *pidx = 0, *sidx = 0, *cidx = 0, *zidx = 0;
  pomp_onestep_sim *ff = NULL;

  if (deltat <= 0)
    errorcall(R_NilValue,"'delta.t' should be a positive number");

    {
      int *dim;
      dim = INTEGER(GET_DIM(xstart)); nvars = dim[0]; nreps = dim[1];
      dim = INTEGER(GET_DIM(params)); npars = dim[0];
      dim = INTEGER(GET_DIM(covar)); covlen = dim[0]; ncovars = dim[1];
      ntimes = LENGTH(times);
    }

    PROTECT(Snames = GET_ROWNAMES(GET_DIMNAMES(xstart))); nprotect++;
    PROTECT(Pnames = GET_ROWNAMES(GET_DIMNAMES(params))); nprotect++;
    PROTECT(Cnames = GET_COLNAMES(GET_DIMNAMES(covar))); nprotect++;

    // set up the covariate table
    struct lookup_table covariate_table = {covlen, ncovars, 0, REAL(tcovar), REAL(covar)};

    // vector for interpolated covariates
    PROTECT(cvec = NEW_NUMERIC(ncovars)); nprotect++;
    SET_NAMES(cvec,Cnames);

    // indices of accumulator variables
    nzeros = LENGTH(zeronames);
    zidx = INTEGER(PROTECT(matchnames(Snames,zeronames,"state variables"))); nprotect++;

    // extract user function
    PROTECT(fn = pomp_fun_handler(func,gnsi,&mode)); nprotect++;

    // set up
    switch (mode) {

    case Rfun:			// R function

      PROTECT(dtvec = NEW_NUMERIC(1)); nprotect++;
      PROTECT(tvec = NEW_NUMERIC(1)); nprotect++;
      PROTECT(xvec = NEW_NUMERIC(nvars)); nprotect++;
      PROTECT(pvec = NEW_NUMERIC(npars)); nprotect++;
      SET_NAMES(xvec,Snames);
      SET_NAMES(pvec,Pnames);

      // set up the function call
      PROTECT(fcall = LCONS(cvec,args)); nprotect++;
      SET_TAG(fcall,install("covars"));
      PROTECT(fcall = LCONS(dtvec,fcall)); nprotect++;
      SET_TAG(fcall,install("delta.t"));
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

    case native:			// native code

      // construct state, parameter, covariate indices
      sidx = INTEGER(PROTECT(matchnames(Snames,GET_SLOT(func,install("statenames")),"state variables"))); nprotect++;
      pidx = INTEGER(PROTECT(matchnames(Pnames,GET_SLOT(func,install("paramnames")),"parameters"))); nprotect++;
      cidx = INTEGER(PROTECT(matchnames(Cnames,GET_SLOT(func,install("covarnames")),"covariates"))); nprotect++;

      *((void **) (&ff)) = R_ExternalPtrAddr(fn);

      break;

    default:

      errorcall(R_NilValue,"unrecognized 'mode' %d",mode); // # nocov

    break;

    }

    // create array to hold results
    {
      int dim[3] = {nvars, nreps, ntimes};
      PROTECT(X = makearray(3,dim)); nprotect++;
      setrownames(X,Snames,3);
    }

    // copy the start values into the result array
    memcpy(REAL(X),REAL(xstart),nvars*nreps*sizeof(double));

    if (mode==1) {
      set_pomp_userdata(args);
      GetRNGstate();
    }

    // now do computations
    {
      int first = 1;
      int use_names = 0;
      int *posn = 0;
      double *time = REAL(times);
      double *xs = REAL(X);
      double *xt = REAL(X)+nvars*nreps;
      double *cp = REAL(cvec);
      double *ps = REAL(params);
      double t = time[0];
      double *pm, *xm;
      int i, j, k, step;

      for (step = 1; step < ntimes; step++, xs = xt, xt += nvars*nreps) {

        R_CheckUserInterrupt();

        if (t > time[step]) {
          errorcall(R_NilValue,"'times' is not an increasing sequence");
        }

        memcpy(xt,xs,nreps*nvars*sizeof(double));

        // set accumulator variables to zero
        for (j = 0; j < nreps; j++)
          for (i = 0; i < nzeros; i++)
            xt[zidx[i]+nvars*j] = 0.0;

        switch (method) {
        case 1:			// one step
          dt = time[step]-t;
          nstep = (dt > 0) ? 1 : 0;
          break;
        case 2:			// fixed step
          dt = deltat;
          nstep = num_map_steps(t,time[step],dt);
          break;
        case 3:			// Euler method
          dt = deltat;
          nstep = num_euler_steps(t,time[step],&dt);
          break;
        default:
          errorcall(R_NilValue,"unrecognized 'method'"); // # nocov
        break;
        }

        for (k = 0; k < nstep; k++) { // loop over Euler steps

          // interpolate the covar functions for the covariates
          table_lookup(&covariate_table,t,cp);

          for (j = 0, pm = ps, xm = xt; j < nreps; j++, pm += npars, xm += nvars) { // loop over replicates

            switch (mode) {

            case Rfun: 		// R function

            {
              double *xp = REAL(xvec);
              double *pp = REAL(pvec);
              double *tp = REAL(tvec);
              double *dtp = REAL(dtvec);
              double *ap;

              *tp = t;
              *dtp = dt;
              memcpy(xp,xm,nvars*sizeof(double));
              memcpy(pp,pm,npars*sizeof(double));

              if (first) {

                PROTECT(ans = eval(fcall,rho));	nprotect++; // evaluate the call
                if (LENGTH(ans) != nvars) {
                  errorcall(R_NilValue,"user 'step.fun' returns a vector of %d state variables but %d are expected: compare initial conditions?",
                            LENGTH(ans),nvars);
                }

                PROTECT(nm = GET_NAMES(ans)); nprotect++;
                use_names = !isNull(nm);
                if (use_names) {
                  posn = INTEGER(PROTECT(matchnames(Snames,nm,"state variables"))); nprotect++;
                }

                ap = REAL(AS_NUMERIC(ans));

                first = 0;

              } else {

                ap = REAL(AS_NUMERIC(PROTECT(eval(fcall,rho))));
                UNPROTECT(1);

              }

              if (use_names) {
                for (i = 0; i < nvars; i++) xm[posn[i]] = ap[i];
              } else {
                for (i = 0; i < nvars; i++) xm[i] = ap[i];
              }

            }

              break;

            case native: 		// native code

              (*ff)(xm,pm,sidx,pidx,cidx,ncovars,cp,t,dt);

              break;

            default:

              errorcall(R_NilValue,"unrecognized 'mode' %d",mode); // # nocov

            break;

            }

          }

          t += dt;

          if ((method == 3) && (k == nstep-2)) { // penultimate step
            dt = time[step]-t;
            t = time[step]-dt;
          }
        }
      }
    }

    if (mode==1) {
      PutRNGstate();
      unset_pomp_userdata();
    }

    UNPROTECT(nprotect);
    return X;
}

int num_euler_steps (double t1, double t2, double *dt) {
  double tol = sqrt(DOUBLE_EPS);
  int nstep;
  // nstep will be the number of Euler steps to take in going from t1 to t2.
  // note also that the stepsize changes.
  // this choice is meant to be conservative
  // (i.e., so that the actual dt does not exceed the specified dt
  // by more than the relative tolerance 'tol')
  // and to counteract roundoff error.
  // It seems to work well, but is not guaranteed:
  // suggestions would be appreciated.

  if (t1 >= t2) {
    *dt = 0.0;
    nstep = 0;
  } else if (t1+*dt >= t2) {
    *dt = t2-t1;
    nstep = 1;
  } else {
    nstep = (int) ceil((t2-t1)/(*dt)/(1+tol));
    *dt = (t2-t1)/((double) nstep);
  }
  return nstep;
}

int num_map_steps (double t1, double t2, double dt) {
  double tol = sqrt(DOUBLE_EPS);
  int nstep;
  // nstep will be the number of discrete-time steps to take in going from t1 to t2.
  nstep = (int) floor((t2-t1)/dt/(1-tol));
  return (nstep > 0) ? nstep : 0;
}
