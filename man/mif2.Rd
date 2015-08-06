\name{Iterated filtering 2}
\title{IF2: Maximum likelihood by iterated, perturbed Bayes maps}
\docType{methods}
\alias{mif2}
\alias{mif2,mif2d.pomp-method}
\alias{mif2-mif2d.pomp}
\alias{mif2,pfilterd.pomp-method}
\alias{mif2-pfilterd.pomp}
\alias{mif2,pomp-method}
\alias{mif2-pomp}
\alias{continue,mif2d.pomp-method}
\alias{continue-mif2d.pomp}
\alias{mif2d.pomp-class}
\alias{mif2d.pomp-methods}
\alias{conv.rec,mif2d.pomp-method}
\alias{conv.rec-mif2d.pomp}
\alias{conv.rec,mif2List-method}
\alias{conv.rec-mif2List}
\alias{coef,mif2List-method}
\alias{coef.rec-mif2List}
\alias{plot-mif2d.pomp}
\alias{plot,mif2d.pomp-method}
\alias{plot-mif2List}
\alias{plot,mif2List-method}
\alias{mif2List-class}
\alias{c-mif2d.pomp}
\alias{c,mif2d.pomp-method}
\alias{c-mif2List}
\alias{c,mif2List-method}
\alias{[-mif2List}
\alias{[,mif2List-method}
\alias{rw.sd}
\description{
  An improved iterated filtering algorithm for estimating the parameters of a partially-observed Markov process.
  Running \code{mif2} causes the algorithm to perform a specified number of particle-filter iterations.
  At each iteration, the particle filter is performed on a perturbed version of the model, in which the parameters to be estimated are subjected to random perturbations at each observation.
  This extra variability effectively smooths the likelihood surface and combats particle depletion by introducing diversity into particle population.
  As the iterations progress, the magnitude of the perturbations is diminished according to a user-specified cooling schedule.
  The algorithm is presented and justified in Ionides et al. (2015).
}
\usage{
\S4method{mif2}{pomp}(object, Nmif = 1, start, Np, rw.sd, transform = FALSE,
    cooling.type = c("hyperbolic", "geometric"), cooling.fraction.50,
    tol = 1e-17, max.fail = Inf, verbose = getOption("verbose"), \dots)
\S4method{mif2}{pfilterd.pomp}(object, Nmif = 1, Np, tol, \dots)
\S4method{mif2}{mif2d.pomp}(object, Nmif, start, Np, rw.sd, transform,
    cooling.type, cooling.fraction.50, tol, \dots)
\S4method{continue}{mif2d.pomp}(object, Nmif = 1, \dots)
\S4method{conv.rec}{mif2d.pomp}(object, pars, transform = FALSE, \dots)
\S4method{conv.rec}{mif2List}(object, \dots)
rw.sd(\dots)
}
\arguments{
  \item{object}{
    An object of class \code{pomp}.
  }
  \item{Nmif}{
    The number of filtering iterations to perform.
  }
  \item{start}{
    named numerical vector;
    the starting guess of the parameters.
    By default, \preformatted{start=coef(object).}
  }
  \item{Np}{
    the number of particles to use in filtering.
    This may be specified as a single positive integer, in which case the same number of particles will be used at each timestep.
    Alternatively, if one wishes the number of particles to vary across timestep, one may specify \code{Np} either as a vector of positive integers (of length \code{length(time(object))}) or as a function taking a positive integer argument.
    In the latter case, \code{Np(n)} must be a single positive integer, representing the number of particles to be used at the \code{n}-th timestep:
    \code{Np(1)} is the number of particles to use going from \code{timezero(object)} to \code{time(object)[1]},
    \code{Np(2)}, from \code{time(object)[1]} to \code{time(object)[2]},
    and so on.
    \strong{Note that this behavior differs from that of \code{\link{mif}}!}
  }
  \item{rw.sd}{
    specification of the magnitude of the random-walk perturbations that will be applied to some or all model parameters.
    Parameters that are to be estimated should have positive perturbations specified here.
    The specification is given using the \code{rw.sd} function, which creates a list of unevaluated expressions.
    The latter are evaluated in a context where the model time variable is defined (as \code{time}).
    The expression \code{ivp(s)} can be used in this context as shorthand for \preformatted{ifelse(time==time[1],s,0).}
    Likewise, \code{ivp(s,lag)} is equivalent to \preformatted{ifelse(time==time[lag],s,0).}
    See below for some examples.
    The perturbations that are applied are normally distributed with the specified s.d.
    If \code{transform = TRUE}, then they are applied on the estimation scale.
  }
  \item{transform}{
    logical;
    if \code{TRUE}, optimization is performed on the estimation scale, as defined by the user-supplied parameter transformations (see \code{\link{pomp}}).
    This can be used, for example, to enforce positivity or interval constraints on model parameters.
    See the tutorials on the \href{http://kingaa.github.io/pomp}{package website} for examples.
  }
  \item{cooling.type, cooling.fraction.50}{
    specifications for the cooling schedule, i.e., the manner in which the intensity of the parameter perturbations is reduced with successive filtering iterations.
    \code{cooling.type} specifies the nature of the cooling schedule.
    See below (under \dQuote{Specifying the perturbations}) for more detail.
  }
  \item{tol, max.fail}{
    passed to the particle filter.
    See the descriptions under \code{\link{pfilter}}.
  }
  \item{verbose}{
    logical; if TRUE, print progress reports.
  }
  \item{\dots}{
    additional arguments that override the defaults.
  }
  \item{pars}{names of parameters.}
}
\value{
  Upon successful completion, \code{mif2} returns an object of class \code{mif2d.pomp}.
  This class inherits from the \code{\link{pfilterd.pomp}} and \code{\link{pomp}} classes.
}
\section{Specifying the perturbations: the \code{rw.sd} function}{
  This function simply returns a list containing its arguments as unevaluated expressions.
  These are then evaluated in a context containing the model \code{time} variable.
  This allows for easy specification of the structure of the perturbations that are to be applied.
  For example,
  \preformatted{
    rw.sd(a=0.05,
          b=ifelse(0.2,time==time[1],0),
          c=ivp(0.2),
          d=ifelse(time==time[13],0.2,0),
          e=ivp(0.2,lag=13),
          f=ifelse(time<23,0.02,0))}
  results in perturbations of parameter \code{a} with s.d. 0.05 at every time step, while parameters \code{b} and \code{c} both get perturbations of s.d. 0.2 only before the first observation.
  Parameters \code{d} and \code{e}, by contrast, get perturbations of s.d. 0.2 only before the thirteenth observation.
  Finally, parameter \code{f} gets a random perturbation of size 0.02 before every observation falling before \eqn{t=23}.

  On the \eqn{m}-th IF2 iteration, prior to time-point \eqn{n}, the \eqn{d}-th parameter is given a random increment normally distributed with mean \eqn{0} and standard deviation \eqn{c_{m,n} \sigma_{d,n}}{c[m,n] sigma[d,n]}, where \eqn{c} is the cooling schedule and \eqn{\sigma}{sigma} is specified using \code{rw.sd}, as described above.
  Let \eqn{N} be the length of the time series and \eqn{\alpha=}{alpha=}\code{cooling.fraction.50}.
  Then, when \code{cooling.type="geometric"}, we have
  \deqn{c_{m,n}=\alpha^{\frac{n-1+(m-1)N}{50N}}.}{c[m,n]=alpha^((n-1+(m-1)N)/(50N)).}
  When \code{cooling.type="hyperbolic"}, we have
  \deqn{c_{m,n}=\frac{s+1}{s+n+(m-1)N},}{c[m,n]=(s+1)/(s+n+(m-1)N),}
  where \eqn{s} satisfies
  \deqn{\frac{s+1}{s+50N}=\alpha.}{(s+1)/(s+50N)=alpha.}
  Thus, in either case, the perturbations at the end of 50 IF2 iterations are a fraction \eqn{\alpha}{alpha} smaller than they are at first.
}
\section{Re-running \code{mif2} Iterations}{
  To re-run a sequence of \code{mif2} iterations, one can use the \code{mif2} method on a \code{mif2d.pomp} object.
  By default, the same parameters used for the original \code{mif2} run are re-used (except for \code{tol}, \code{max.fail}, and \code{verbose}, the defaults of which are shown above).
  If one does specify additional arguments, these will override the defaults.
}
\section{Continuing \code{mif2} Iterations}{
  One can resume a series of \code{mif2} iterations from where one left off using the \code{continue} method.
  A call to \code{mif2} to perform \code{Nmif=m} iterations followed by a call to \code{continue} to perform \code{Nmif=n} iterations will produce precisely the same effect as a single call to \code{mif2} to perform \code{Nmif=m+n} iterations.
  By default, all the algorithmic parameters are the same as used in the original call to \code{mif2}.
  Additional arguments will override these defaults.
}
\section{Methods}{
  Methods that can be used to manipulate, display, or extract information from a \code{mif2d.pomp} object:
  \describe{
    \item{conv.rec}{
      \code{conv.rec(object, pars = NULL)} returns the columns of the convergence-record matrix corresponding to the names in \code{pars}.
      By default, all rows are returned.
    }
    \item{logLik}{
      Returns the value in the \code{loglik} slot.
      NB: this is \emph{not} the same as the likelihood of the model at the MLE!
    }
    \item{c}{
      Concatenates \code{mif2d.pomp} objects into a \code{mif2List}.
    }
    \item{plot}{
      Plots a series of diagnostic plots when applied to a \code{mif2d.pomp} or \code{mif2List} object.
    }
  }
}
\examples{
\dontrun{
pompExample(ou2)

guess1 <- guess2 <- coef(ou2)
guess1[c('x1.0','x2.0','alpha.2','alpha.3')] <- 0.5*guess1[c('x1.0','x2.0','alpha.2','alpha.3')]
guess2[c('x1.0','x2.0','alpha.2','alpha.3')] <- 1.5*guess1[c('x1.0','x2.0','alpha.2','alpha.3')]

m1 <- mif2(ou2,Nmif=100,start=guess1,Np=1000,
           cooling.type="hyperbolic",cooling.fraction.50=0.05,
           rw.sd=rw.sd(x1.0=ivp(0.5),x2.0=ivp(0.5),
             alpha.2=0.1,alpha.3=0.1))

m2 <- mif2(ou2,Nmif=100,start=guess2,Np=1000,
           cooling.type="hyperbolic",cooling.fraction.50=0.05,
           rw.sd=rw.sd(x1.0=ivp(0.5),x2.0=ivp(0.5),
             alpha.2=0.1,alpha.3=0.1))

plot(c(m1,m2))

rbind(mle1=c(coef(m1),loglik=logLik(pfilter(m1,Np=1000))),
      mle2=c(coef(m2),loglik=logLik(pfilter(m1,Np=1000))),
      truth=c(coef(ou2),loglik=logLik(pfilter(m1,Np=1000))))
}
}
\references{
  E. L. Ionides, D. Nguyen, Y. Atchad{\\'e}, S. Stoev, and A. A. King. 
  Inference for dynamic and latent variable models via iterated, perturbed Bayes maps.
  Proc. Natl. Acad. Sci. U.S.A., 112:719--724, 2015.
}
\author{Aaron A. King, Edward L. Ionides, and Dao Nguyen}
\seealso{
  \code{\link{pomp}}, \code{\link{pfilter}}, \code{\link{mif}}, and the tutorials on the \href{http://kingaa.github.io/pomp}{package website}.
}
\keyword{optimize}
\keyword{ts}
