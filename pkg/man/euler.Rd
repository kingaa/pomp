\name{euler}
\alias{euler}
\alias{euler.simulate}
\alias{euler.density}
\title{Dynamical models based on stochastic Euler algorithms}
\description{
  Facilities for simulating discrete-time Markov processes and continuous-time Markov processes using the Euler algorithm.
}
\usage{
euler.simulate(xstart, times, params, step.fun, delta.t, \dots,
               statenames = character(0), paramnames = character(0),
               covarnames = character(0), zeronames = character(0),
               tcovar, covar, PACKAGE)
euler.density(x, times, params, dens.fun, \dots,
              statenames = character(0), paramnames = character(0),
              covarnames = character(0), tcovar, covar, log = FALSE,
              PACKAGE)
}
\arguments{
  \item{xstart}{
    Matrix (dimensions \code{nvar} x \code{nrep}) of states at initial time \code{times[1]}.
  }
  \item{x}{
    Matrix (dimensions \code{nvar} x \code{nrep} x \code{ntimes}) of states at times \code{times}.
  }
  \item{times}{
    Vector of times (length \code{ntimes}) at which states are required or given.
  }
  \item{params}{
    Matrix containing parameters of the model.
    The \code{nrep} columns of \code{params} correspond to those of \code{xstart}.
  }
  \item{step.fun}{
    This can be either an R function or a compiled, dynamically loaded native function containing the model simulator.
    It should be written to take a single Euler step from a single point in state space.
    If it is a native function, it must be of type "euler\_step\_sim" as defined in the header "pomp.h", which is included with the package.
    For details on how to write such codes, see Details.
  }
  \item{dens.fun}{
    This can be either an R function or a compiled, dynamically loaded native function containing the model transition probability density function.
    This function will be called to compute the probability of the actual Euler steps.
    It must be of type "euler\_step\_pdf" as defined in the header "pomp.h", which is included with the package.
    For details on how to write such codes, see Details.
  }
  \item{delta.t}{
    Time interval of Euler steps.
  }
  \item{statenames}{
    Names of state variables, in the order they will be expected by the routine named in \code{euler.step.fun} and \code{euler.dens.fun}.
  }
  \item{paramnames}{
    Names of parameters, in the order they will be expected by the routine named in \code{euler.step.fun} and \code{euler.dens.fun}.
  }
  \item{covarnames}{
    Names of columns of the matrix of covariates \code{covar}, in the order they will be expected by the routine named in \code{euler.step.fun} and \code{euler.dens.fun}.
  }
  \item{zeronames}{
    Names of additional variables which will be zeroed before each time in \code{times}.
    These are useful, e.g., for storing accumulations of state variables.
  }
  \item{tcovar}{
    Times at which covariates are measured.
  }
  \item{covar}{
    Matrix of covariates.
    This should have dimensions \code{length(tcovar)} x \code{ncovar}, where \code{ncovar} is the number of covariates.
  }
  \item{log}{
    logical; if TRUE, probabilities p are given as log(p).
  }
  \item{\dots}{
    if \code{step.fun} (or \code{dens.fun}) is an R function, then additional arguments will be passed to it.
    If \code{step.fun} (or \code{dens.fun}) is a native routine, then additional arguments are ignored.
  }
  \item{PACKAGE}{
    an optional argument that specifies to which dynamically loaded library we restrict the search for the native routines.
    If this is '"base"', we search in the R executable itself.
  }
}
\details{
  If \code{step.fun} is written as an R function, it must have at least the arguments \code{x}, \code{t}, \code{params}, \code{delta.t}, and \code{\dots}.
  On a call to this function, \code{x} will be a named vector of state variables, \code{t} a scalar time, and \code{params} a named vector of parameters.
  The length of the Euler step will be \code{delta.t}.
  If the argument \code{covars} is included and a covariate table has been included in the \code{pomp} object, then on a call to this function, \code{covars} will be filled with the values, at time \code{t}, of the covariates.
  This is accomplished via interpolation of the covariate table.
  Additional arguments may be given: these will be filled by the correspondingly-named elements in the \code{userdata} slot of the \code{pomp} object (see \code{\link{pomp}}).

  If \code{step.fun} is written in a native language, it must be a function of type "euler\_step\_sim" as specified in the header "pomp.h" included with the package (see the directory "include" in the installed package directory).

  If \code{dens.fun} is written as an R function, it must have at least the arguments \code{x1}, \code{x2}, \code{t1}, \code{t2}, \code{params}, and \code{\dots}.
  On a call to this function, \code{x1} and \code{x2} will be named vectors of state variables at times \code{t1} and \code{t2}, respectively.
  The named vector \code{params} contains the parameters.
  If the argument \code{covars} is included and a covariate table has been included in the \code{pomp} object, then on a call to this function, \code{covars} will be filled with the values, at time \code{t1}, of the covariates.
  This is accomplished via interpolation of the covariate table.
  As above, any additional arguments will be filled by the correspondingly-named elements in the \code{userdata} slot of the \code{pomp} object (see \code{\link{pomp}}).

  If \code{dens.fun} is written in a native language, it must be a function of type "euler\_step\_pdf" as defined in the header "pomp.h" included with the package (see the directory "include" in the installed package directory).
}
\value{
  \code{euler.simulate} returns a \code{nvar} x \code{nrep} x \code{ntimes} array, where \code{nvar} is the number of state variables, \code{nrep} is the number of replicate simulations (= number of columns of \code{xstart} and \code{params}), and \code{ntimes} is the length of \code{times}.
  If \code{x} is this array, \code{x[,,1]} will be identical to \code{xstart}; the rownames of \code{x} and \code{xstart} will also coincide.

  \code{euler.density} returns a \code{nrep} x \code{ntimes-1} array.
  If \code{f} is this array, \code{f[i,j]} is the likelihood of a transition from \code{x[,i,j]} to \code{x[,i,j+1]} in exactly one Euler step of duration \code{times[j+1]-times[j]}.
}
\author{Aaron A. King (kingaa at umich dot edu)}
\seealso{\code{\link{eulermultinom}}, \code{\link{pomp}}}
\examples{
## set up a lookup table for basis functions for the seasonality
tbasis <- seq(0,25,by=1/52)
basis <- periodic.bspline.basis(tbasis,nbasis=3)
colnames(basis) <- paste("seas",1:3,sep='')

## some parameters
params <- c(
            gamma=26,mu=0.02,iota=0.01,
            beta1=1200,beta2=1800,beta3=600,
            beta.sd=1e-3,
            pop=2.1e6,
            rho=0.6,
            S.0=26/1200,I.0=0.001,R.0=1-0.001-26/1200
            )

## set up the pomp object
## the C codes "sir_euler_simulator", "sir_euler_density", and "sir_ODE" are included
## in the "examples" directory (file "sir.c")
po <- pomp(
           time=seq(1/52,20,by=1/52),
           data=rbind(measles=numeric(52*20)),
           t0=0,
           tcovar=tbasis,
           covar=basis,
           delta.t=1/52/20,
           statenames=c("S","I","R","cases","W","B","dW"),
           paramnames=c("gamma","mu","iota","beta1","beta.sd","pop","rho"),
           covarnames=c("seas1"),
           zeronames=c("cases"),
           step.fun="sir_euler_simulator",
           dens.fun="sir_euler_density",
           skeleton="sir_ODE",
           rprocess=euler.simulate,
           dprocess=euler.density,
           dmeasure="binom_dmeasure",
           rmeasure="binom_rmeasure",
           initializer=function(params,t0,...){
             p <- exp(params)
             with(
                  as.list(p),
                  {
                    fracs <- c(S.0,I.0,R.0)
                    x0 <- c(
                            round(pop*fracs/sum(fracs)), # make sure the three compartments sum to 'pop' initially
                            rep(0,9)	# zeros for 'cases', 'W', and the transition numbers
                            )
                    names(x0) <- c("S","I","R","cases","W","B","SI","SD","IR","ID","RD","dW")
                    x0
                  }
                  )
           }
           )

## simulate from the model
tic <- Sys.time()
x <- simulate(po,params=log(params),nsim=3)
toc <- Sys.time()
print(toc-tic)
plot(x[[1]])
}
\keyword{models}
