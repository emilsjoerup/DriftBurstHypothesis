library(testthat)
library(DriftBurstHypothesis)

context("DBH C++ test")
test_that("DBH error codes",{
  #these functions has broken before, hopefully they won't break again.
  expect_equal(DriftBurstHypothesis:::AsymptoticVarianceC(c(1:3), 3), NaN) 
  expect_equal(DriftBurstHypothesis:::AsymptoticVarianceC(c(1:3), 4), NaN)
  expect_equal(DriftBurstHypothesis:::AutomaticLagSelectionC(1:10, 30) , 7)
  }
)

context("DBH test")
test_that("DBH sim test", {
  set.seed(123)
  iT = 23400; Mean_bandwidth = 300L
  timestamps = seq(0, 23400, length.out = iT+1)
  testtimes  = seq(0, 23400, 60L)
  
  r = rnorm(iT, mean = 0.02, sd = 1)/sqrt(iT)
  p = c(0,cumsum(r))
  
  
  DBH = drift_bursts(timestamps, p, testtimes, PreAverage = 1, AcLag = -1, 
                     Mean_bandwidth = Mean_bandwidth, Variance_bandwidth = 5*Mean_bandwidth, bParallelize = FALSE)
  expect_equal(mean(getMu(DBH)[-1]),0.02516096)
  expect_equal(mean(getSigma(DBH)[-1]), 0.9629604)
})


context("Examples check")
test_that("DBH Examples check",{
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
  testtimes = seq(34200, 57600, 60)
  logprices = log(vY)
  
  library("DriftBurstHypothesis")
  
  #calculating drift burst hypothesis
  
  DBH = drift_bursts(timestamps,  logprices,
                     testtimes, PreAverage = 5, AcLag = -1L,
                     Mean_bandwidth = 300L, Variance_bandwidth = 900L,
                     bParallelize = FALSE)
  
  
  #plot test statistic
  plot = plot(DBH)
  #plot both test statistic and price
  plot2 = plot(DBH, price = vY, time = timestamps)
  #Plot the mu series
  plot3 = plot(DBH, which = "Mu")
  #plot the sigma series
  plot4 = plot(DBH, which = "Sigma")
  expect_equal(plot3$usr[1:2], plot4$usr[1:2])
  #plot both
  plot5 = plot(DBH, which = c("Mu", "Sigma"))
  
  #Means of the tstat, sigma, and mu series.
  expect_equal(mean(getDB(DBH)), 0.01859701)
  expect_equal(mean(getSigma(DBH)), 0.12483728)
  expect_equal(mean(getMu(DBH)), - 0.001151098)
  
  
  
  
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
  StartTime = strptime("1970-01-01 00:00:00.0000", "%Y-%m-%d %H:%M:%OS", tz = "GMT")
  Tradetime = StartTime + timestamps
  testtimes = seq(34200, 57600, 60)
  
  
  price = xts(vY, Tradetime)
  
  
  DBHxts = drift_bursts(time = NULL,  log(price), 
                        testtimes, PreAverage = 5, AcLag = -1L,
                        Mean_bandwidth = 300L, Variance_bandwidth = 900L, 
                        bParallelize = FALSE)
  plot6 = plot(DBHxts)
  expect_true(all.equal.list(plot, plot6))
  plot7 = plot(DBHxts, price = price)
  expect_true(all.equal.list(plot2, plot7))
  #check for equality
  expect_true(all.equal(as.numeric(getDB(DBH)), as.numeric(getDB(DBHxts))))
})