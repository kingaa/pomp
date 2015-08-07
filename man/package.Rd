\name{pomp-package}
\title{Inference for partially observed Markov processes}
\docType{package}
\alias{The pomp package}
\alias{pomp-package}
\alias{pomp package}
\description{
  The \pkg{pomp} package provides facilities for inference on time series data using partially-observed Markov process (\acronym{POMP}) models.
  These models are also known as state-space models or nonlinear stochastic dynamical systems.
  One can use \pkg{pomp} to fit nonlinear, non-Gaussian dynamic models to time-series data.
  The package is both a set of tools for data analysis and a platform upon which statistical inference methods for \acronym{POMP} models can be implemented.
}
\section{Data analysis using \pkg{pomp}}{
  The first step in using \pkg{pomp} is to encode one's model(s) and data in objects of class \code{pomp}.
  One does this via a call to \code{\link{pomp}}, which involves specifying the unobserved state process and the measurement process of the model.
  Details on this are given in the documentation for the \code{\link{pomp}} constructor function.
  Examples are given in the tutorials on the \href{http://kingaa.github.io/pomp}{package website}, in the demos (\code{demo(package='pomp')}), and via the \code{\link{pompExample}} function.

  \pkg{pomp} version \Sexpr[echo=FALSE,stage=build,results=text]{packageDescription("pomp",fields="Version")} provides algorithms for
  \enumerate{
    \item simulation of stochastic dynamical systems; see \code{\link[=simulate-pomp]{simulate}}
    \item particle filtering (AKA sequential Monte Carlo or sequential importance sampling); see \code{\link{pfilter}}
    \item the iterated filtering methods of Ionides et al. (2006, 2011, 2015);
    see \code{\link{mif}} and \code{\link{mif2}}
    \item the nonlinear forecasting algorithm of Kendall et al. (2005); see \code{\link{nlf}}
    \item the particle MCMC approach of Andrieu et al. (2010); see \code{\link{pmcmc}}
    \item the probe-matching method of Kendall et al. (1999, 2005); see \code{\link{probe.match}}
    \item a spectral probe-matching method (Reuman et al. 2006, 2008); see \code{\link{spect.match}}
    \item synthetic likelihood a la Wood (2010); see \code{\link{probe}}
    \item approximate Bayesian computation (Toni et al. 2009); see \code{\link{abc}}
    \item the approximate Bayesian sequential Monte Carlo scheme of Liu & West (2001); see \code{\link{bsmc}}
    \item simple trajectory matching; see \code{\link{traj.match}}.
  }
  The package also provides various tools for plotting and extracting information on models and data.
}
\section{Developing inference tools on the \pkg{pomp} platform}{
  \pkg{pomp} provides a very general interface to the components of \acronym{POMP} models.
  All the inference algorithms in \pkg{pomp} interact with the models and data via this interface.
  One goal of the \pkg{pomp} project has been to facilitate the development of new algorithms in an environment where they can be tested and compared on a growing body of models and datasets.
  
  The low-level interface relevant to developers is documented \link[=rprocess]{here}.
}
\section{Comments, bug reports, feature requests}{
  Contributions are welcome, as are comments, feature requests, and bug reports.
  See the package website \url{http://kingaa.github.io/pomp} for more information, access to the package new RSS feed, links to the authors' websites, references to the literature, and up-to-date versions of the package source and documentation.
}
\section{Documentation}{
  A number of tutorials, demonstrating the construction of \code{pomp} objects and the application of various inference algorithms, are available on the package homepage: \url{http://kingaa.github.io/pomp}.
  Several examples of the construction of \code{pomp} objects are provided with the package and are documented in the help pages: to view a full list of these, execute \code{pompExample()}.
  In addition, there are a number of demos, which can be viewed by executing \code{demo(package="pomp")}.
}    
\section{History}{
  Much of the groundwork for \pkg{pomp} was laid by a working group of the National Center for Ecological Analysis and Synthesis (\acronym{NCEAS}), \dQuote{Inference for Mechanistic Models}.
}
\section{License}{
  \pkg{pomp} is provided under the \acronym{GNU} Public License (\acronym{GPL}).
}
\references{
  See the package website, \url{http://kingaa.github.io/pomp}, for the references.
}
\author{Aaron A. King}
\seealso{
  \code{\link{pomp}},
  \link{pomp low-level interface}, 
  \code{\link{pfilter}},
  \code{\link[=simulate-pomp]{simulate}},
  \code{\link{mif}},
  \code{\link{nlf}},
  \code{\link{probe}},
  \code{\link{traj.match}},
  \code{\link{bsmc2}},
  \code{\link{pmcmc}}
}
\keyword{models}
\keyword{datasets}
\keyword{ts}
