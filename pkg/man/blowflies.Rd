\name{blowflies}
\alias{blowflies}
\alias{blowflies1}
\alias{blowflies2}
\docType{data}
\title{Model for Nicholson's blowflies.}
\description{
  \code{blowflies1} and \code{blowflies2} are \code{pomp} objects encoding stochastic delay-difference models.
}
\usage{data(blowflies)}
\examples{
data(blowflies)
plot(blowflies1)
plot(blowflies2)
}
\details{
  The data are from "population I", a control culture in one of A. J. Nicholson's experiments with the Australian sheep-blowfly \emph{Lucilia cuprina}.
  The experiment is described on pp. 163--4 of Nicholson (1957).
  Unlimited quantities of larval food were provided;
  the adult food supply (ground liver) was constant at 0.4g per day.
  The data were taken from the table provided by Brillinger et al. (1980).
}
\seealso{\code{\link{pomp-class}} and the vignettes}
\references{
  A. J. Nicholson (1957)
  The self-adjustment of populations to change.
  Cold Spring Harbor Symposia on Quantitative Biology, \bold{22}, 153--173.
  
  Y. Xia and H. Tong (2011)
  Feature Matching in Time Series Modeling.
  \emph{Statistical Science} \bold{26}, 21--46.

  E. L. Ionides (2011)
  Discussion of ``Feature Matching in Time Series Modeling'' by Y. Xia and H. Tong.
  \emph{Statistical Science} \bold{26}, 49--52.

  S. N. Wood (2010)
  Statistical inference for noisy nonlinear ecological dynamic systems.
  \emph{Nature} \bold{466}, 1102--1104.
 
  W. S. C. Gurney, S. P. Blythe, and R. M. Nisbet (1980)
  Nicholson's blowflies revisited.
  \emph{Nature} \bold{287}, 17--21.

  D. R. Brillinger, J. Guckenheimer, P. Guttorp and G. Oster (1980)
  Empirical modelling of population time series: The case of age and density dependent rates.
  in G. Oster (ed.), Some Questions in Mathematical Biology, vol. 13, pp. 65--90.
  American Mathematical Society, Providence.
}
\keyword{datasets}
