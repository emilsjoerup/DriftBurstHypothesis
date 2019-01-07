\name{drift_bursts}
\alias{drift_bursts}
\title{
Drift Bursts
}
\description{
Calculates the Test-Statistic for the Drift Burst Hypothesis.
}
\usage{
drift_bursts(time, logprices, testtimes = seq(34200, 57600, 60),
             PreAverage = 5, AcLag = -1L, Mean_bandwidth = 300L, 
             Variance_bandwidth = 900L, bParallelize = FALSE, iCores = NA)
}

\arguments{
  \item{time}{
A \code{numeric} of timestamps for the trades. The timestamps used for testing were of the format second after midnight of the given day, but the precision was in nano-seconds. 
}
  \item{logprices}{
A \code{numeric} containing the logarithm of the prices of the asset.
}
  \item{testtimes}{
A \code{numeric} containing the times at which to calculate the tests. The standard of \code{seq(floor(min(time)), ceiling(max(time)), 60)} denotes calculating the test-statistic once per minute, i.e. 391 times for a typical 6.5 hour trading day.
}
  \item{PreAverage}{
An \code{integer} denoting the period for pre-averaging estimates of the log-prices.
}
  \item{AcLag}{
An \code{integer} denoting which lag is to be used for the HAC estimator of the variance - the standard of \code{-1} denotes using an automatic lag selection algorithm.
}
  \item{Mean_bandwidth}{
An \code{integer} denoting the bandwidth for the exponential kernel for the mean.
}
  \item{Variance_bandwidth}{
An \code{integer} denoting the bandwidth for the exponential kernel for the variance.
}
  \item{bParallelize}{
A \code{Boolean} to determine whether to parallelize the underlying C++ code. 
}
\item{iCores}{
An \code{integer} denoting the number of cores to use for calculating the code when paralellized.
  }
}
\details{
The DBH test statistic cannot be calculated before 1 period of testtimes has passed.

The test statistic is unstable before \code{max(Mean_bandwidth , Variance_bandwidth)}
}
\value{
A list containing the series of the drift burst hypothesis test-statistic as well as the estimated mean and variance series.
}
\references{
Christensen, Oomen and Reno (2018) <DOI:10.2139/ssrn.2842535>.
}
\author{
Emil Sjoerup
}

\examples{
#Simulate from a stochastic volatility model.
#A flash crash  coded into the function.
StochasticVolatilitySim = function(iT, dSigma, dPhi, dMu){
  vSim = numeric(iT)
  vEps = rnorm(iT , sd =dSigma)
  vEpsy = rnorm(iT)
  vEps[30001:32000] = rnorm(2000 ,sd =seq(from = dSigma , 
                                          to = 2*dSigma , length.out = 2000)) 
  vEps[32001:34000] = rnorm(2000 ,sd =seq(from = 2*dSigma , 
                                          to = dSigma , length.out = 2000))
  vEpsy[30001:32000] = -rnorm(2000 ,mean =seq(from = 0,
                                              to = 0.4 , length.out = 2000)) 
  vEpsy[32001:34000] = -rnorm(2000 ,mean =seq(from = 0.4,
                                              to = 0 , length.out = 2000))

  vSim[1] = dMu + dPhi *rnorm(1 , mean = dMu , sd = dSigma /sqrt(1-dPhi^2))
  for (i in 2:iT) {
    vSim[i] = dMu + dPhi * (vSim[(i-1)] - dMu) + vEps[i]
  }
  vY = exp(vSim/2) * vEpsy
  return(vY)
}

iT = 46000; dSigma = 0.4; dPhi = 0.98; dMu = -10;

set.seed(123)
vY = 500+cumsum(StochasticVolatilitySim(iT, dSigma, dPhi, dMu))

#insert an outlier to illustrate robustness.
vY[45000] = 500

#Here, the observations are equidistant, but the code can handle 
#unevenly spaced observations.
timestamps = seq(34200 , 57600 , length.out = iT)
testtimes = seq(34200, 57600, 60)
logprices = log(vY)
DBHtStat = drift_bursts(timestamps,  logprices, 
                        testtimes, PreAverage = 5, AcLag = -1L,
                        Mean_bandwidth = 300L, Variance_bandwidth = 900L, 
                        bParallelize = FALSE)

#plot series
plot(vY , x = timestamps/86400 , type = 'l')
#plot test statistic
plot(DBHtStat$DriftBursts, x = testtimes/86400, type = 'l')

}