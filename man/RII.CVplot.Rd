\name{RII.CVplot}
\alias{RII.CVplot}
\title{RII Cross Validation Plots}
\description{
Evaluates and plots a cross validation score over a specified range of values of the smoothing parameter.
}
\usage{
RII.CVplot(count, pop, loglambda, ...)
}
\arguments{
  \item{count}{ A matrix of outcome counts with number of rows equal to the number of classes and number of columns equal to the number of standardizing groups.  The outcome might be death, or disease incidence, for example. }
  \item{pop}{ A matrix of amounts of exposure, with dimension the same as that of \code{count}.  The amount of exposure could be, for example, the number of person-years at risk, the mid-study period population or the number of individuals at risk at the start of the study period. }
  \item{loglambda}{ A vector of values of the smoothing parameter (on log scale) }
  \item{...}{ Graphical parameters }
}
\details{  
\code{RII.CVplot} is designed to be used for exploratory purposes before estimating a relative index of inequality using \code{RII}.  The cross validation score is computed taking each element of \code{loglambda} in turn to be the value of the smoothing parameter and the results are plotted.  Using \code{RII.CVplot} (with perhaps several different values of \code{loglambda}) should allow the user to choose a suitable argument \code{grid} for use in \code{RII}, or select a single value of the smoothing parameter that is close to being optimal.
}
\references{ Sergeant, J. C. and Firth D. (2004)  Relative index of inequality: definition, estimation and inference.  In preparation. }
\author{Jamie Sergeant, \email{jamie.sergeant@nuffield.oxford.ac.uk} }
\seealso{\code{\link{RII}}, \code{\link{plot.RII}}.}
\examples{
## Plot the cross validation score over a range
## of smoothing parameter values for the LSDeaths data
data(LSDeaths)
LSdead <- xtabs(Deaths ~ class + age, data = LSDeaths)
LSatrisk <- xtabs(AtRisk ~ class + age, data = LSDeaths)
RII.CVplot(LSdead, LSatrisk, loglambda = seq(-2,18,len=21))
}
\keyword{models}
\keyword{regression}
