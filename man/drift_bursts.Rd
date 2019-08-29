\name{drift_bursts}
\alias{drift_bursts}

\title{
deprecated - currently a wrapper for driftBursts
}
\description{
This is a deprecated wrapper for the \link{driftBursts} function.
}
\usage{
drift_bursts(time = NULL, logprices, testtimes = seq(34200, 57600, 60),
                        PreAverage = 5, AcLag = -1L, Mean_bandwidth = 300L, 
                        Variance_bandwidth = 900L, bParallelize = FALSE, iCores = NA,
                        warnings = TRUE)
}
\arguments{
  \item{time}{
See \link{driftBursts}.
}
  \item{logprices}{
See \link{driftBursts}.
}
  \item{testtimes}{
See \link{driftBursts}.
}
  \item{PreAverage}{
See \link{driftBursts}.

}
  \item{AcLag}{
See \link{driftBursts}.
}
  \item{Mean_bandwidth}{
See \link{driftBursts}.
}
  \item{Variance_bandwidth}{
See \link{driftBursts}.
}
  \item{bParallelize}{
See \link{driftBursts}.
}
\item{iCores}{
See \link{driftBursts}.
}
\item{warnings}{
See \link{driftBursts}.
}
}
\details{

See \link{driftBursts}.

}
\value{
See \link{driftBursts}.
}
\references{
Christensen, Oomen and Reno (2018) <DOI:10.2139/ssrn.2842535>.
}
\author{
Emil Sjoerup
}

\examples{
#This function is deprecated, please see help(driftBursts).
}