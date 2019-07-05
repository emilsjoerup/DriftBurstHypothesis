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
  expect_equal( mean(getMu(DBH)[-1]),0.02516096)
  expect_equal(mean(getSigma(DBH)[-1]), 0.9629604)
})
