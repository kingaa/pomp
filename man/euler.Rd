\name{euler}
\alias{euler}
\alias{euler.simulate}
\alias{euler.density}
\title{Dynamical models based on stochastic Euler algorithms}
\description{
  Facilities for simulating stochastic Markov processes based on the Euler algorithm.
}
\usage{
euler.simulate(xstart, times, params, euler.step.fun, delta.t,
               statenames, paramnames, zeronames, tcovar, covar,
               PACKAGE)
euler.density(x, times, params, euler.dens.fun,
              statenames = character(0), paramnames = character(0),
              tcovar, covar, log = FALSE, PACKAGE)
}
\arguments{
  \item{xstart}{
    Matrix (dimensions \code{nvar} x \code{nsim}) of states at initial time \code{times[1]}.
  }
  \item{x}{
    Matrix (dimensions \code{nvar} x \code{nsim} x \code{ntimes}) of states at times \code{times}.
  }
  \item{times}{
    Vector of times (length \code{ntimes}) at which states are required or given.
  }
  \item{params}{
    Matrix containing parameters of the model.
    The \code{nsim} columns of \code{params} correspond to those of \code{xstart}.
  }
  \item{euler.step.fun}{
    Name of a compiled, loaded native function containing the model simulator.
    This function will be called to take the actual Euler steps.
    It must be of type "euler\_step\_sim" as defined in the header "euler.h", which is included with the package.
  }
  \item{euler.dens.fun}{
    Name of a compiled, loaded native function containing the model probability density function.
    This function will be called to compute the probability of the actual Euler steps.
    It must be of type "euler\_step\_pdf" as defined in the header "euler.h", which is included with the package.
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
  \item{PACKAGE}{
    an optional argument that specifies to which dynamically loaded library we restrict the search for \code{euler.step.fun}.
    If this is '"base"', we search in the R executable itself.
  }
}
\value{
  \code{euler.simulate} returns a \code{nvar} x \code{nsim} x \code{ntimes} array, where \code{nvar} is the number of state variables, \code{nsim} is the number of replicate simulations (= number of columns of \code{xstart} and \code{params}), and \code{ntimes} is the length of \code{times}.
  If \code{x} is this array, \code{x[,,1]} will be identical to \code{xstart}; the rownames of \code{x} and \code{xstart} will also coincide.

  \code{euler.density} returns a \code{nsim} x \code{ntimes-1} array.
  If \code{f} is this array, \code{f[i,j]} is the likelihood of a transition from \code{x[,i,j]} to \code{x[,i,j+1]} in exactly one Euler step of duration \code{times[j+1]-times[j]}.
}
\author{Aaron A. King (kingaa at umich dot edu)}
\seealso{\code{\link{eulermultinom}}, \code{\link{pomp}}}
\examples{
## set up a lookup table for basis functions for the seasonality
tbasis <- seq(0,20,by=1/52)
basis <- periodic.bspline.basis(tbasis,nbasis=3)

## some parameters
params <- c(
            gamma=26,mu=0.2,iota=0.01,
            beta1=1200,beta2=2100,beta3=300,
            beta.sd=0.1,
            pop=2.1e5,
            rho=0.6,
            S.0=26/1200,I.0=0.001,R.0=1-0.001-26/1200
            )

## set up the pomp object
## the C codes "sir_euler_simulator" and "sir_euler_density" are included in the "examples" directory
po <- pomp(
           time=seq(1/52,20,by=1/52),
           data=rbind(measles=numeric(52*20)),
           t0=0,
           tbasis=tbasis,
           basis=basis,
           dt=1/52/20,
           statenames=c("S","I","R","cases","W","trans1"),
           paramnames=c("gamma","mu","iota","beta1","beta.sd","pop"),
           zeronames=c("cases"),
           rprocess=function(xstart,times,params,dt,tbasis,basis,
             statenames,paramnames,zeronames,...){
             euler.simulate(
                            xstart=xstart,
                            times=times,
                            params=params,
                            euler.step.fun="sir_euler_simulator",
                            delta.t=dt,
                            statenames=statenames,
                            paramnames=paramnames,
                            zeronames=zeronames,
                            tcovar=tbasis,
                            covar=basis
                            )
           },
           dprocess=function(x,times,params,log,tbasis,basis,
             statenames,paramnames,...){
             euler.density(
                           x=x,
                           times=times,
                           params=params,
                           euler.dens.fun="sir_euler_density",
                           statenames=statenames,
                           paramnames=paramnames,
                           tcovar=tbasis,
                           covar=basis,
                           log=log
                           )
           },
           rmeasure=function(x,times,params,...){
             nsims <- ncol(params)
             ntimes <- length(times)
             array(
                   data=rbinom(
                     n=nsims*ntimes,
                     size=x["cases",,],
                     prob=exp(params["rho",])
                     ),
                   dim=c(1,nsims,ntimes),
                   dimnames=list("measles",NULL,NULL)
                   )
           },
           dmeasure=function(y,x,times,params,log=FALSE,...){
             nreps <- ncol(params)
             ntimes <- length(times)
             f <- array(dim=c(nreps,ntimes))
             for (k in 1:nreps)
               f[k,] <- dbinom(
                               x=y,
                               size=x["cases",k,],
                               prob=exp(params["rho",k]),
                               log=log
                               )
             f
           },
           initializer=function(params,t0,...){
             p <- exp(params)
             pop <- p["pop"]
             fracs <- p[c("S.0","I.0","R.0")]
             x0 <- c(
                     round(pop*fracs/sum(fracs)), # make sure the three compartments sum to 'pop' initially
                     rep(0,9)	# zeros for 'cases', 'W', and the transition numbers
                     )
             names(x0) <- c(
                            "S","I","R","cases","W",
                            paste("trans",1:7,sep="")
                            )
             x0
           }
           )

## simulate from the model
tic <- Sys.time()
x <- simulate(po,params=log(params),nsim=3)
toc <- Sys.time()
print(toc-tic)

t <- seq(0,4/52,1/52/20)
X <- simulate(po,params=log(params),nsim=10,states=TRUE,obs=TRUE,times=t)

f <- dprocess(
              po,
              x=X$states[,,31:40],
              times=t[31:40],
              params=matrix(
                log(params),
                nrow=length(params),
                ncol=10,
                dimnames=list(names(params),NULL)
                ),
              log=TRUE
              )
apply(f,1,sum)

g <- dmeasure(
              po,
              y=t(as.matrix(X$obs[,7,])),
              x=X$states,
              times=t,
              params=matrix(
                log(params),
                nrow=length(params),
                ncol=10,
                dimnames=list(names(params),NULL)
                ),
              log=TRUE
              )
apply(g,1,sum)
}
\keyword{models}
