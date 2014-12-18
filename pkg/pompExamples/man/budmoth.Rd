\name{budmoth}
\docType{data}
\alias{budmoth.sim}
\title{Larch budmoth model POMPs with real and simulated data.}
\description{
  \code{pompExample(budmoth.sim)} constructs a \code{pomp} object containing the larch budmoth model and simulated budmoth density, parasitism rate, and food quality (needle-length) data.
  Four datasets, representing four distinct parameter regimes, are available.

  The model has three state variables: 
  \enumerate{
    \item{\eqn{Q_t} (measure of food quality on \eqn{[0,1]})} 
    \item{\eqn{N_t} (budmoth density)}
    \item{\eqn{S_t} (fraction of budmoth larvae infected with parasitoids)}
  }

  There are three observables:
  \enumerate{
    \item{\eqn{\hat Q_t} (needle length)}
    \item{\eqn{\hat N_t}}
    \item{\eqn{\hat S_t}}
  }
}
\section{State process}{
  Uncorrelated random effects, for \eqn{t=1,\dots,T}:
  \deqn{\alpha_t \sim \mathrm{LogitNormal}(\mathrm{logit}(\alpha),\sigma_{\alpha}^2)}
  \deqn{\lambda_t \sim \mathrm{Gamma}(\lambda,\sigma_{\lambda}^2)}
  \deqn{a_t \sim \mathrm{LogNormal}(\log(a),\sigma_{a}^2)}
  Note: \eqn{X} is \eqn{\mathrm{LogitNormal}(\mu,\sigma)} if \eqn{\mathrm{logit}(X)} is \eqn{\mathrm{Normal}(\mu,\sigma)}. 
  
  The inverse of \eqn{\mathrm{logit}} is \eqn{\mathrm{expit}}.
  \pkg{R} functions \code{logit}, \code{expit}, \code{rlogitnorm}, \code{dlogitnorm} are defined in \pkg{pompExamples}.
  
  The state process, for \eqn{t=1,\dots,T}:
  \deqn{Q_{t} = (1-\alpha_{t})\frac{\gamma}{\gamma+N_{t-1}} +\alpha_{t}Q_{t-1}}
  \deqn{N_{t} = \lambda_t N_{t-1} (1-S_{t-1})\exp\big\{-gN_{t-1}-\delta(1-Q_{t-1})\big\}}
  \deqn{S_{t} =  1-\exp\left(\frac{-a_tS_{t-1}N_{t-1}}{1+a_twS_{t-1}N_{t-1}} \right)}
}  
\section{Measurement process}{
  For \eqn{t=1,\dots,T}:
  \deqn{\hat Q_t \sim \mathrm{LogNormal}(\log(\beta_0+\beta_1Q_t),\sigma_Q)}
  \deqn{\hat N_t \sim \mathrm{LogNormal}(\log(N_t),\sigma_N)}
  \deqn{\hat S_t \sim \mathrm{LogitNormal}(\mathrm{logit}(uS_t),\sigma_S)}
}
\examples{
## three regimes, high and low noise regimes for parasitism and tritrophic
\dontrun{
bm <- pompExample(budmoth,envir=NULL)
plot(bm$food)
plot(bm$para1)
plot(bm$para2)
plot(bm$tri)
}
}
\keyword{models}
