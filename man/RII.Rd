\name{RII}
\alias{RII}
\alias{print.RII}
\alias{summary.RII}
\alias{print.summary.RII}
\title{Relative Index of Inequality Estimation}
\description{
Estimates the relative index of inequality (Sergeant and Firth, 2004) for data consisting of the amount of exposure and observed numbers of outcomes in a set of ordered (socio-economic) classes, possibly cross-classified by some standardizing variable such as age.  The (continuous) incidence rate is modelled by fitting smoothing splines by maximum penalized likelihood, with smoothing parameter selection by cross validation.
}
\usage{
RII(count, pop, loglambda = NULL, grid = NULL, se = FALSE, B = NULL, alpha = 0.025, returnboot = FALSE)
}
\arguments{
  \item{count}{ A matrix of outcome counts with number of rows equal to the number of classes and number of columns equal to the number of standardizing groups.  The outcome might be death, or disease incidence, for example. }
  \item{pop}{ A matrix of amounts of exposure, with dimension the same as that of \code{count}.  The amount of exposure could be, for example, the number of person-years at risk, the mid-study period population or the number of individuals at risk at the start of the study period. }
  \item{loglambda}{ Optional value of the smoothing parameter (on log scale) }
  \item{grid}{ A vector of values (on log scale) on which to search for a starting value for use in the optimization of the smoothing parameter }
  \item{se}{ Should a bootstrap standard error be computed? }
  \item{B}{ The number of bootstrap samples to use if computing a standard error }
  \item{alpha} { If a standard error is computed, an approximate 1-2*\code{alpha} confidence interval for the RII is also computed}
  \item{returnboot} { If a standard error is computed, should the bootstrap datasets and their respective RII estimates and values of \code{loglambda} be returned?}
}
\details{  
  If \code{loglambda} is supplied then this value of the smoothing parameter is used in \bold{all} calculations and \code{grid} is redundant.  For no smoothing, specify \code{loglambda = -Inf}.  Specifying \code{loglambda = Inf} will induce a linear fit.

If \code{loglambda} is not supplied then \code{grid} is required.  The element of \code{grid} which yields the smallest value of the cross validation score is taken to be the starting value in a minimization of the score over the smoothing parameter.  If this element is equal to \code{min(grid)}, the optimum \code{loglambda} is taken to be \code{-Inf}.  If this element is equal to \code{max(grid)}, the optimum \code{loglambda} is taken to be \code{Inf}.  If bootstrapping is performed, \code{grid} is used for each bootstrap dataset.

For a given value of the smoothing parameter, the penalized Poisson log likelihood is maximized.
}
\value{
An object of class \code{RII}, with some of the components
  \item{count}{\code{count}}
  \item{pop}{\code{pop}}
  \item{loglambda}{The value of the smoothing parameter (on log scale) used to estimate the RII}
  \item{par}{The optimum spline coefficients.  When \code{loglambda = Inf} these are the intercept and gradient of the linear fit.}
  \item{group.effects}{Standardizing group effects}
  \item{maxval}{Maximum value of the penalized log likelihood}
  \item{expected}{The fitted outcome counts}
  \item{residuals}{Deviance residuals}
  \item{var}{Delta method approximation of var(\code{RII})}
  \item{var.log}{Delta method approximation of var(log(\code{RII}))}
  \item{RII}{The estimated RII}
  \item{se}{Bootstrap standard error for log(\code{RII})}
  \item{alpha}{\code{alpha}}
  \item{ci}{Approximate 1-2*\code{alpha} empirical percentile interval for the RII}
  \item{boot.data}{The \code{B} bootstrap datasets}
  \item{boot.rep}{The estimated RIIs for the bootstrap datasets}
  \item{boot.lambda}{The values of \code{loglambda} used to estimate the RIIs for the bootstrap datasets}
}
\references{ Sergeant, J. C. and Firth D. (2004)  Relative index of inequality: definition, estimation and inference.  In preparation. }
\note{
Methods available for objects of class \code{RII} are 
\itemize{
  \item \code{plot.RII}
  \item \code{print.RII}
  \item \code{summary.RII}
  \item \code{print.summary.RII} }
}
\author{Jamie Sergeant, \email{jamie.sergeant@nuffield.oxford.ac.uk} }
\seealso{\code{\link{plot.RII}}, \code{\link{RII.CVplot}}.}
\examples{
## Estimate the RII for the LSDeaths data,
## using a smoothing parameter of 1
data(LSDeaths)
LSdead <- xtabs(Deaths ~ class + age, data = LSDeaths)
LSatrisk <- xtabs(AtRisk ~ class + age, data = LSDeaths)
LSRII <- RII(LSdead, LSatrisk, loglambda = 0)
}
\keyword{models}
\keyword{regression}
