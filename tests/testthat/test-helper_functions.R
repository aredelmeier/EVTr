library(testthat)
library(EVTr)

context("test-helper_functions.R")

# For using same Random number generator as CircelCI (R version 3.5.x)
RNGversion(vstr = "3.5.0")

test_that("Test functions in explanation.R", {

  data <- data.frame(ID = c(rep(1, 3), rep(2, 2), rep(3, 4), 4), obs_id = c(1, 2, 3, 1, 2, 1, 2, 3, 4, 1), Injury_Length = rexp(10))
  data <- data.table(data)
  exc <- excess(data, ne = 1)

  # must include data as a data.frame
  expect_error(excess(data = 1, ne = 1))

  # data must include column called "Injury_data"
  expect_error(excess(data = data.frame(ID = c(1, 2, 3)), ne = 1))

  # must include "ne" - number above the threshold
  expect_error(excess(data))

  # error when all data below threshold
  expect_error(excess(data, ne = 100))

  # error when all data above threshold
  expect_error(excess(data, ne = -100))

  # expect less rows of data when looking above threshold
  expect_gt(nrow(data), nrow(exc))

  # threshold must be less than original data
  for (i in 1:nrow(exc)) {
    expect_gte(exc$Injury_Length_before[i], exc$threshold[1])
  }

  # Injury_Length = Injury_Length_before - threshold
  for (i in 1:nrow(exc)) {
    expect_equal(exc$Injury_Length[i], exc$Injury_Length_before[i] - exc$threshold[1])
  }

})
