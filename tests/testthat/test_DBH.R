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


