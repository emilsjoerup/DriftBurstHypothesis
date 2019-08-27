# DriftBurstHypothesis

An R-package for the calculation of the Drift Burst Hypothesis test-statistic from the working paper Christensen, Oomen and Reno (2018) <DOI:10.2139/ssrn.2842535>.

The t-statistic at period n is calculated as follows:

![equation](https://latex.codecogs.com/png.latex?%5Cbar%7BT%7D%5En%20%3D%20%5Csqrt%7B%5Cfrac%7Bh_%7Bn%7D%7D%7BK_%7B2%7D%7D%7D%5Cfrac%7B%5Chat%7B%5Cbar%7B%5Cmu%7D%7D_%7Bt%7D%5E%7Bn%7D%7D%7B%5Csqrt%7B%5Chat%7B%5Cbar%7B%5Csigma%7D%7D_%7Bt%7D%5E%7Bn%7D%7D%7D), 

where the local mean estimator is:

![equation](https://latex.codecogs.com/png.latex?%5Chat%7B%5Cbar%7B%5Cmu%7D%7D_%7Bt%7D%5E%7Bn%7D%3D%5Cfrac%7B1%7D%7Bh_%7Bn%7D%7D%5Csum_%7Bi%3D1%7D%5E%7Bn-k_%7Bn%7D&plus;2%7DK%5Cleft%28%5Cfrac%7Bt_%7Bi-1%7D-t%7D%7Bh_%7Bn%7D%7D%5Cright%29%5CDelta_%7Bi-1%7D%5E%7Bn%7D%5Coverline%7BY%7D),

and the local variance estimator is:

![equation](https://latex.codecogs.com/png.latex?%5Chat%7B%5Cbar%7B%5Csigma%7D%7D_%7Bt%7D%5E%7Bn%7D%20%3D%20%5Cfrac%7B1%7D%7Bh_%7Bn%7D%27%7D%5Cleft%5B%5Csum_%7Bi%3D1%7D%5E%7Bn-k_%7Bn%7D&plus;2%7D%5Cleft%28K%5Cleft%28%5Cfrac%7Bt_%7Bi-1%7D-t%7D%7Bh%27_%7Bn%7D%7D%5Cright%29%5CDelta_%7Bi-1%7D%5E%7Bn%7D%5Coverline%7BY%7D%5Cright%29%5E%7B2%7D&plus;2%5Csum_%7BL%3D1%7D%5E%7BL_%7Bn%7D%7D%5Comega%5Cleft%28%5Cfrac%7BL%7D%7BL_%7Bn%7D%7D%5Cright%29%5Csum_%7Bi%3D1%7D%5E%7Bn-k_%7Bn%7D-L&plus;2%7DK%5Cleft%28%5Cfrac%7Bt_%7Bi-1%7D-t%7D%7Bh_%7Bn%7D%27%7D%5Cright%29K%5Cleft%28%5Cfrac%7Bt_%7Bi&plus;L-1%7D-t%7D%7Bh_%7Bn%7D%27%7D%5Cright%29%5CDelta_%7Bi-1%7D%5E%7Bn%7D%5Coverline%7BY%7D%5CDelta_%7Bi-1&plus;L%7D%5E%7Bn%7D%5Coverline%7BY%7D%5Cright%5D)


with:

![equation](https://latex.codecogs.com/png.latex?%5CDelta_%7Bi%7D%5E%7Bn%7D%5Coverline%7BY%7D%20%3D%20%5Csum_%7Bj%3D1%7D%5E%7Bk_%7Bn%7D-1%7Dg_%7Bj%7D%5E%7Bn%7D%5CDelta_%7Bi&plus;j%7D%5E%7Bn%7DY)

denoting the overlapping pre-averaged returns with the weighting function:


![equation](https://latex.codecogs.com/png.latex?g%5Cleft%28x%5Cright%29%3D%5Ctext%7Bmin%7D%5Cleft%28x%2C1-x%5Cright%29),

and

![equation](https://latex.codecogs.com/png.latex?%5Comega%5Cleft%28%5Ccdot%5Cright%29)

is a smooth kernel defined on the positive real numbers, ![equation](https://latex.codecogs.com/png.latex?L_%7Bn%7D) is the lag length over which the estimator is applied. By default, the lag-length will be determined by way of the Newey-West algorithm.





## Examples using simulated high frequency data:
```
library(highfrequency)
library(xts)
library(DriftBurstHypothesis)
data("sample_tdata")
price = xts(as.numeric(sample_tdata$PRICE), index(sample_tdata))
plot(price)
testtimes = seq(34200, 57600, 60)



DBHxts = drift_bursts(timestamps = NULL,  logpricexts,
                      testTimes, preAverage = 5, ACLag = -1L,
                      meanBandwidth = 300L, varianceBandwidth = 900L,
                      bParallelize = TRUE, iCores = 8)

plot(DBHxts, price = price)
```

![Example plot](https://i.imgur.com/LTWi7uM.png)


```
library(DriftBurstHypothesis)
set.seed(1234)
returns = rnorm(23399, sd = 1)/sqrt(23400)
price = c(0,cumsum(returns))
timestamps = seq(34200, 57600, length.out = 23400)
testTimes = c(34200,seq(34200 + 5*300, 57600, 60))
DBH = driftBursts(timestamps, price, testTimes, preAverage = 5, meanBandwidth = 300, varianceBandwidth = 5*300)
plot(DBH, price = price, timestamps = timestamps)
```

![Example plot](https://imgur.com/a/RI0S2Pi)