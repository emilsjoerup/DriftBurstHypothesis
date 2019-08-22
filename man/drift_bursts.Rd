\name{drift_bursts}
\alias{drift_bursts}
\alias{getMu}
\alias{getSigma}
\alias{getDB}
\title{
Drift Bursts
}
\description{
Calculates the Test-Statistic for the Drift Burst Hypothesis.
}
\usage{
drift_bursts(timestamps = NULL, logPrices, testTimes = seq(34200, 57600, 60),
                        preAverage = 5, ACLag = -1L, meanBandwidth = 300L, 
                        varianceBandwidth = 900L, bParallelize = FALSE,
                        iCores = NA, warnings = TRUE)
}

\arguments{
  \item{timestamps}{
Either:
A \code{numeric} of timestamps for the trades in seconds after midnight.
Or: NULL, when the time argument is NULL, the logprices argument must be an \code{xts} object.
}
  \item{logPrices}{
A \code{numeric} or \code{xts} object containing the log-prices of the asset.
}
  \item{testTimes}{
A \code{numeric} containing the times at which to calculate the tests. The standard of \code{seq(34200, 57600, 60)} denotes calculating the test-statistic once per minute, i.e. 391 times for a typical 6.5 hour trading day from 9:30:00 to 16:00:00.
See details.
}
  \item{preAverage}{
An \code{integer} denoting the length of pre-averaging estimates of the log-prices.

}
  \item{ACLag}{
An positive \code{integer} denoting how many lags are to be used for the HAC estimator of the variance - the standard of \code{-1} denotes using an automatic lag selection algorithm for each iteration..
}
  \item{meanBandwidth}{
An \code{integer} denoting the bandwidth for the left-sided exponential kernel for the mean.
}
  \item{varianceBandwidth}{
An \code{integer} denoting the bandwidth for the left-sided exponential kernel for the variance.
}
  \item{bParallelize}{
A \code{Boolean} to determine whether to parallelize the underlying C++ code. (Using OpenMP)
}
\item{iCores}{
An \code{integer} denoting the number of cores to use for calculating the code when paralellized.
If this argument is not provided, sequential evaluation will be used even though \code{bParallelize} is TRUE
}
\item{warnings}{
A \code{logical} denoting whether warnings should be shown.
}
}
\details{

If the testtimes vector contains instructions to test before the first trade (excluding the very first entry), or more than 15 minutes after the last trade, these entries will be deleted, as not doing so may cause crashes.

The test statistic is NOT calculated on the first entry of testtimes, please use this to denote the start of a trading day.

The test statistic is unstable before \code{max(meanBandwidth , varianceBandwidth)} seconds has passed.

If times is provided and logprices is an \code{xts} object, the indices of logprice will be used regardless.

Note that using an \code{xts} argument is slower than using a \code{numeric} and the timestamps due to the creation of the timestamps from the index of the input.


The lags from the Newey-West algorithm is increased by 2 * (PreAveage-1) as due to the Pre-Averaging we know atleast this many lags should be corrected for due to the auto-correlation. The maximum of 20 lags is also increased by this same factor for the same reason.


}
\value{
An object of class \code{DBH} and \code{list} containing the series of the drift burst hypothesis test-statistic as well as the estimated spot drift and variance series. The list also contains some information such as the variance and mean bandwidths along with the pre-averaging setting and the amount of observations.

The S3 class DBH has the following methods so far:
\itemize{
\item{\code{getDB}: Extracts the test-statistic.}
\item{\code{getSigma}: Extracts the estimated local variance. Note that for preAverage != 1, this value cannot be directly interpreted as the volatility as it can with preAverage = 1, a correction for this is being worked on. This does not affect the test statistic!}
\item{\code{getMu}: Extracts the estimated local drift. Note that for preAverage != 1, this value cannot be directly interpreted as the drift as it can with preAverage = 1, a correction for this is being worked on. This does not affect the test statistic!}
\item{\code{plot}: Plotting routine.}
The plotting routine accepts 10 arguments for customization in the ... argument, which must be passed as a named list:
\itemize{
  \item{\code{which}:Used for choosing which series to plot, valid choices are: \code{"DriftBursts"}, \code{"DB"}, \code{"Sigma"}, \code{"Mu"}, and     \code{c("Sigma","Mu")}, the order of the latter is irrelevant. The case of the input does not matter. Default = \code{"DriftBursts"}}
  \item{\code{price}:The price series which, if provided, will be overlayed in a red dotted line and the level will be shown at the right y-axis.       (Only   used if which is \code{"DriftBursts"}). Default = \code{NULL}}
  \item{\code{time}:Timestamps for the trades in seconds after midnight, which will be used for the x-axis of the plot if the price is overlayed.       Default = \code{NULL}}
  \item{\code{startTime}:Start of the trading day in seconds after midnight. Default = \code{34200}}
  \item{\code{endTime}:End of the trading day in seconds after midnight. Default = \code{57600}}
  \item{\code{xlab}:The label of the x-axis. Default = \code{"time"}}
  \item{\code{ylab}:The label of the y-axis. Default = \code{"value"}}
  \item{\code{leg.x}:X-position of the legend in the case \code{which} is \code{"DriftBursts"} AND the price is overlayed. Default = \code{"topleft"}.   Usage is as in the base \code{R} engine.}
  \item{\code{leg.y}:Y-position of the legend in the case \code{which} is \code{"DriftBursts"} AND the price is overlayed. Default = \code{"NULL"}.     Usage is as in the base \code{R} engine.}
  \item{\code{tz}: \code{character} denoting the time-zone. Default = \code{"GMT"}}
  \item{\code{annualize}: \code{logical} denoting whether to annualize in case which contains "mu", "sigma", or both. Default = FALSE}
  \item{\code{nDays}: \code{numeric} denoting how many to use for annualization if \code{annualize} is TRUE. Default = 252}
  }
}
}
\references{
Christensen, Oomen and Reno (2018) <DOI:10.2139/ssrn.2842535>.
}
\author{
Emil Sjoerup
}

\examples{
#Simulate from a stochastic volatility model.
#Both a flash crash and flash rally are coded into the function.
StochasticVolatilitySim = function(iT, dSigma, dPhi, dMu){
  vSim = numeric(iT)
  vEps = rnorm(iT , sd =dSigma)
  vEpsy = rnorm(iT)
  vEps[30001:32000] = rnorm(2000 ,sd =seq(from = dSigma , 
                                          to = 2*dSigma , length.out = 2000)) 
  vEps[32001:34000] = rnorm(2000 ,sd =seq(from = 2*dSigma , 
                                          to = dSigma , length.out = 2000))
  vEpsy[30001:32000] = -rnorm(2000 ,mean =seq(from = 0,
                                              to = 0.3 , length.out = 2000)) 
  vEpsy[32001:34000] = -rnorm(2000 ,mean =seq(from = 0.3,
                                              to = 0 , length.out = 2000))
  
  
  vEps[60001:63000] = rnorm(3000,sd = seq(from = dSigma , 
                                          to = 2*dSigma , length.out = 3000))
  vEps[63001:66000] = rnorm(3000,  sd = seq(from = 2*dSigma , 
                                            to =  dSigma, length.out = 3000))
  
  vEpsy[60001:63000] = rnorm(3000 ,mean =seq(from = 0,
                                             to = 0.2 , length.out = 3000))
  vEpsy[63001:66000] = rnorm(3000 ,mean =seq(from = 0.2,
                                             to = 0 , length.out = 3000))
  vSim[1] = dMu + dPhi *rnorm(1 , mean = dMu , sd = dSigma /sqrt(1-dPhi^2))
  for (i in 2:iT) {
    vSim[i] = dMu + dPhi * (vSim[(i-1)] - dMu) + vEps[i]
  }
  vY = exp(vSim/2) * vEpsy
  return(vY)
}
#Set parameter values of the simulation
iT = 66500; dSigma = 0.3; dPhi = 0.98; dMu = -10;
#set seed for reproducibility
set.seed(123)
#Simulate the series
vY = 500+cumsum(StochasticVolatilitySim(iT, dSigma, dPhi, dMu))

#insert an outlier to illustrate robustness.
vY[50000] = 500

#Here, the observations are equidistant, but the code can handle unevenly spaced observations.
timestamps = seq(34200 , 57600 , length.out = iT)
testTimes = seq(34200, 57600, 60)
logprices = log(vY)

library("DriftBurstHypothesis")

#calculating drift burst hypothesis

DBH = drift_bursts(timestamps,  logprices,
                        testTimes, preAverage = 5, ACLag = -1L,
                        meanBandwidth = 300L, varianceBandwidth = 900L,
                        bParallelize = FALSE)


#plot test statistic
plot(DBH)
#plot both test statistic and price
plot(DBH, price = vY, time = timestamps)
#Plot the mu series
plot(DBH, which = "Mu")
#plot the sigma series
plot(DBH, which = "Sigma")
#plot both
plot(DBH, which = c("Mu", "Sigma"))


#Means of the tstat, sigma, and mu series.
mean(getDB(DBH))
mean(getSigma(DBH))
mean(getMu(DBH))



\dontrun{
################## same example with xts object:
library("xts")
#Set parameter values of the simulation
iT = 66500; dSigma = 0.3; dPhi = 0.98; dMu = -10;
#set seed for reproducibility
set.seed(123)
#Simulate the series
vY = 500+cumsum(StochasticVolatilitySim(iT, dSigma, dPhi, dMu))

#insert an outlier to illustrate robustness.
vY[50000] = 500

#Here, the observations are equidistant, but the code can handle unevenly spaced observations.
timestamps = seq(34200 , 57600 , length.out = iT)
StartTime = strptime("1970-01-01 00:00:00.0000", "\%Y-\%m-\%d \%H:\%M:\%OS", tz = "GMT")
Tradetime = StartTime + timestamps
testTimes = seq(34200, 57600, 60)


price = xts(vY, Tradetime)


DBHxts = drift_bursts(timestamps = NULL,  log(price), 
                   testTimes, preAverage = 5, ACLag = -1L,
                   meanBandwidth = 300L, varianceBandwidth = 900L, 
                   bParallelize = FALSE)

#check for equality
all.equal(as.numeric(getDB(DBH)), as.numeric(getDB(DBHxts)))
}
}