library(testthat)
library(EVTr)

context("test-helper_functions.R")

# For using same Random number generator as CircelCI (R version 3.5.x)
# RNGversion(vstr = "3.5.0")

test_that("Test functions in gen_data.R", {
  censor <- 5
  xi <- 0.5
  n <- 10
  num_inj <- 3
  rate_exp <- 1
  ne <- 5

  data1 <- gen_data(censor = censor, xi = xi, n = n, num_inj = num_inj, rate_exp = rate_exp, specific =
                      "delete_censored_obs")

  # No observations should be censored
  expect_true(all(data1$Censored == FALSE))

  # simulates n individuals
  expect_true(max(data1$ID) == n)

  # maximum num_inj number of rows per ID
  data1[, num := .N, by = ID]
  expect_true(max(data1$num) <= num_inj)

  data2 <- gen_data(censor = censor, xi = xi, n = n, num_inj = num_inj, rate_exp = rate_exp, specific =
                      "keep_censored_obs")

  # more rows then the dt that removes the cesnored observations
  expect_gte(nrow(data2), nrow(data1))

  # calculate proportion of censored injuries correctly
  expect_equal(data2$Prop_of_injuries_censored[1], sum(data2$Censored) / n)

  # calculate proportion of censored injuries correctly
  expect_equal(data2$Prop_of_ind_with_censored[1], sum(data2$Censored) / nrow(data2))

  # all censored injuries must be less than the original injury lengths
  expect_true(all(data2$Actual >= data2$Injury_Length))

  data3 <- gen_data(censor = censor, xi = xi, n = n, num_inj = num_inj, rate_exp = rate_exp, specific =
                      "keep_only_max_obs")

  # we extract one observation for each individual
  expect_true(nrow(data3) == n)

  # we calculate the proportion correctly
  expect_equal(sum(data3$Censored) / n, data3$Prop_maxima_censored[1])


  data4 <- gen_data(censor = censor, xi = xi, n = n, num_inj = num_inj, rate_exp = rate_exp, ne = ne, specific =
                      "max_excess")

  # we extract at most one observation for each individual
  expect_true(nrow(data4) <= n)

  # we calculate the proportion correctly
  expect_equal(nrow(data4) / n, data4$Prop_maxima_above_thresh[1])

  # Original injuries longer than excesses
  expect_true(all(data4$Actual >= data4$Injury_Length))

  data5 <- gen_data(censor = censor, xi = xi, n = n, num_inj = num_inj, rate_exp = rate_exp, ne = ne, specific =
                      "excess")

  # more rows when you don't take the maxima
  expect_true(nrow(data5) >= nrow(data4))

  # Original injuries longer than excesses
  expect_true(all(data5$Actual >= data5$Injury_Length))

})
