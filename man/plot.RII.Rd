\name{plot.RII}
\alias{plot.RII}
\title{Plot RII Objects}
\description{
A plot method for objects of class \code{RII}.
}
\usage{
\method{plot}{RII}(x, group = 1, ylim = NULL, ...)
}
\arguments{
  \item{x}{ An object of class \code{RII}}
  \item{group}{ Which standardizing group to plot }
  \item{ylim}{ Graphical parameter with a useful default }
  \item{...}{ Graphical parameters }
}
\details{
\code{RII.plot} gives a graphical interpretation of the output from \code{RII}.  For the specified standardizing group, the empirical incidence rate for each class is obtained by dividing the count in that class by the amount of eposure and then plotted.  The natural cubic spline incidence rate fitted by \code{RII} (multiplied by the exponential of the appropriate group effect) is then added to the plot.  For the first standardizing group (and for data with only one group) the group effect is defined to be 0.
}
\references{ Sergeant, J. C. and Firth D. (2004)  Relative index of inequality: definition, estimation and inference.  In preparation. }
\author{Jamie Sergeant, \email{jamie.sergeant@nuffield.oxford.ac.uk} }
\seealso{\code{\link{RII}}, \code{\link{RII.CVplot}}.}
\examples{
## Estimate the RII for the LSDeaths data,
## using a smoothing parameter of 1
data(LSDeaths)
LSdead <- xtabs(Deaths ~ class + age, data = LSDeaths)
LSatrisk <- xtabs(AtRisk ~ class + age, data = LSDeaths)
LSRII <- RII(LSdead, LSatrisk, loglambda = 0)

## Plot the data and fitted death rate for each age group
par(mfrow=c(2,2))
plot(LSRII, group = "25-34", main = "(a) age group 25-34")
plot(LSRII, group = "35-44", main = "(b) age group 35-44")
plot(LSRII, group = "45-54", main = "(c) age group 45-54")
plot(LSRII, group = "55-64", main = "(d) age group 55-64")
}
\keyword{models}
\keyword{regression}
