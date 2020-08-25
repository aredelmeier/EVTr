library(testthat)
library(EVTr)

context("test-estimation.R")

# For using same Random number generator as CircelCI (R version 3.5.x)
RNGversion(vstr = "3.5.0")

test_that("Test mle function in estimation.R", {

  set.seed(1)
  data <- data.frame(Injury_Length = rexp(10), Censored = rbinom(10, 1, 0.5))

  set.seed(1)
  data2 <- data.frame(Injury_Length = rexp(10), Censored = rep(0, 10))

  mle(data = data, method = "MLE")

  # must include proper data - otherwise won't work with optim
  expect_error(mle(data = data.frame(1), method = "MLE"))

  # must include proper MLE method otherwise will throw an error
  expect_error(mle(data = data, method = "MLE222"))

  # the same data should give the same estimates
  expect_equal(mle(data = data, method = "MLE"), mle(data = data2, method = "MLE"))

  # the same data but different mle method should give different estimates
  expect_true(mle(data = data, method = "MLE") != mle(data = data, method = "CensMLE"))

})
