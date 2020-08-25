
#' Simulate injuries according to a generalized Pareto distribution and healthy periods according to a exponential distribution.
#' Based on a fixed parameter \code{censor}, some of the injuries and healthy periods will be censored based on the
#' cummulative sum of the healthy and injury periods. Based on \code{specific}, certain observations will be returned.
#'
#' @details
#' \code{delete_censored_obs}: simulates data, censors the observations, but only returns the fully observed (injury)
#' periods.
#' \code{keep_censored_obs}: simulates data, censors the observations, and returns all (injury) periods.
#' \code{keep_only_max_obs}: simulates data, censors, returns only the longest injury (censored or not censored)
#' for each individual.
#' \code{max_excess}: simulates data, censors, calculates the longest injury (censored or not censored), then
#' calculates the threshold based on \code{ne} (number of exceedances above the threshold). Finally returns only
#' the excess above this given threshold.
#' \code{excess}: simulates data, censors, then calculates the threshold based on \code{ne} (number of exceedances above
#' the threshold). Finally returns only the excess above this given threshold.
#'
#' @param censor Integer. Where to censor observations.
#'
#' @param xi Numeric. The xi parameter in the generalized Pareto distribution.
#'
#' @param n Integer. The number of individuals.
#'
#' @param num_inj Integer. The number of injuries per individual.
#'
#' @param rate_exp The rate in the exponential distribution.
#'
#' @param ne Integer. Number of observations above the threshold. Only used when \code{specific} = 'max_excess' or
#' 'excess'.
#'
#' @param specific String. Corresponding to one of "delete_censored_obs", "keep_censored_obs", "keep_only_max_obs",
#' "max_excess" or "excess".
#'
#' @return Dataframe.
#'
#' This dataframe will have various columns depending on \code{specific}.
#'
#' If \code{specific} = "delete_censored_obs", the first column in \code{data_return} is "Injury_Length" and corresponds
#' to the length of the injury.
#'
#' The second column "ID" stands for the invididuals ID or identifier. "Censored" = 1 if the
#' injury is censored before the entire injury is finished (in which case "Injury_Length" is the length of the injury until
#' the censoring) and = 0 otherwise.
#'
#' "Actual" indicates the length of the injury before it was censored. In other words, if there was no censoring, we would
#' see the actual length of the injury.
#'
#' "Any_Injury_Censored" == 1 if the individual has had their last injury censored. Since all individuals are censored at one
#' point, each individual while either have an injury censored or a healthy period censored. But since we don't care about
#' modelling healthy periods, if the healthy period is censored, we will have full information on the injuries.
#'
#' "Prop_of_injuries_censored" is the # injuries censored / # individuals (n).
#'
#' "Prop_of_ind_with_censored_injury" is the # injuries censored / total # injuries.
#'
#' "Prop_max_censored" is the # individuals where the censored injury is their longest injury / # individuals.
#'
#' If \code{specific} = "keep_censored_obs", the columns are the same.
#'
#' If \code{specific} = "keep_only_max_obs", the above columns all apply. A new column "delta" is:
#' 1 - (1 - Censored) * Any_Injury_Censored
#'
#' If \code{specific} = "max_excess", the above columns all apply.  In this dataframe, "Injury_Length" is the excess injury length
#' above the threshold decided by \code{ne}. A new column "Injury_Length_before" is the original length of the injury (i.e
#' threshold - Injury_Length).
#'
#' "Prop_of_obs_above_threshold" = # injuries >= threshold / total # of injuries.
#'
#' If \code{specific} = "excess", the columns are the same as "max_excess".
#'
#' @export
#'
#' @author Annabelle Redelmeier
#'
#' @import data.table
#'
#' @examples
#' censor <- 3; xi <- 0.5; n <- 10; num_inj <- 3; rate_exp <- 1; ne <- 5
#'
#' gen_data(censor = censor, xi = xi, n = n, num_inj = num_inj, rate_exp = rate_exp, specific = "delete_censored_obs")[]
#'
#' gen_data(censor = censor, xi = xi, n = n, num_inj = num_inj, rate_exp = rate_exp, specific = "keep_censored_obs")[]
#'
#' gen_data(censor = censor, xi = xi, n = n, num_inj = num_inj, rate_exp = rate_exp, specific = "keep_only_max_obs")[]
#'
#' gen_data(censor = censor, xi = xi, n = n, num_inj = num_inj, rate_exp = rate_exp, ne = ne, specific = "max_excess")[]
#'
#' gen_data(censor = censor, xi = xi, n = n, num_inj = num_inj, rate_exp = rate_exp, ne = ne, specific = "excess")[]
#'

gen_data <- function (censor = NULL,
                      xi,
                      n,
                      num_inj,
                      rate_exp,
                      ne = NULL,
                      specific = c("delete_censored_obs",
                                   "keep_censored_obs",
                                   "keep_only_max_obs",
                                   "max_excess",
                                   "excess"),
                      seed = 1) {
  set.seed(seed)

  # Generate data
  Healthy <- matrix(rexp(n * num_inj, rate = rate_exp), nrow = n, ncol = num_inj)
  Injured <- matrix(QRM::rGPD(n*num_inj, xi = xi), nrow = n, ncol = num_inj)

  # Combine matrices
  df <- matrix(NA, nrow = n, ncol = 2 * num_inj)
  df[, seq(from = 1, to = num_inj * 2, by = 2)] <- Healthy
  df[, seq(from = 2, to = num_inj * 2, by = 2)] <- Injured

  df <- as.data.table(df)
  df[, ID := .I]

  df_long <- melt(df, id.vars = "ID", value.name = "Injury_Length")
  setkeyv(df_long, "ID")
  df_long$csum <- ave(df_long$Injury_Length, df_long$ID, FUN = cumsum)

  df_long[, obs_id := 1:.N, by = ID]
  df_long[, variable := NULL]

  df_long[, injury := ifelse(obs_id %% 2 == 0, 1, 0)]


  if (!is.null(censor)) {
    df_long[, Censored := ifelse(csum <= censor, 0, 1)]

    df_long[, lagged_censored :=  shift(Censored), by = ID]

    df_long[, first_censored := ifelse((Censored == 1 & lagged_censored == 0), 1, 0)]
    df_long[, after_censored := ifelse((Censored == 1 & lagged_censored == 1), 1, 0)]

    df_rm <- df_long[after_censored == 0]
    df_rm[, lagged_censored := NULL]
    df_rm[, after_censored := NULL]

    #
    df_rm[, Censored_Length := ifelse(csum > censor, csum - censor, 0)]
    df_rm[, Injury_Length_with_censored := ifelse(Censored == 0, Injury_Length, Injury_Length - Censored_Length)]

    df_rm[, Actual := Injury_Length]
    df_rm[, Injury_Length := Injury_Length_with_censored]
    df_rm[, Injury_Length_with_censored := NULL]
    df_rm[, csum := NULL]
    df_rm[, Censored_Length := NULL]

    data <- df_rm


  } else {
    df_long[, Censored := 0]
    df_long[, first_censored := 0]
    df_long[, Actual := Injury_Length]
    df_long[, Any_Injury_Censored := 0]
    data <- df_long

  }

  Prop_injuries_censored <- nrow(data[injury == 1][first_censored == 1]) / n
  Prop_of_ind_with_censored_injury <- nrow(data[injury == 1][first_censored == 1]) / nrow(data[injury == 1])
  # Prop_healthy_censored <- nrow(data[injury == 0][first_censored == 1]) / n
  # Prop_censored <- nrow(data[first_censored == 1]) / n

  data_inj <- data[injury == 1][, max_injury_length := max(Injury_Length), by = ID]
  tmp <- data_inj[max_injury_length == Injury_Length][Censored == 1]
  Prop_max_censored <- nrow(tmp) / n


  data[, 'Prop_of_injuries_censored' := Prop_injuries_censored]
  data[, 'Prop_of_ind_with_censored_injury' := Prop_of_ind_with_censored_injury]
  data[, 'Prop_maxima_censored' := Prop_max_censored]

  #data[, 'Prop_of_healthy_censored' := Prop_healthy_censored]
  #data[, 'Prop_of_injuries_and_healthy_censored' := Prop_censored]

  # Extract only injuries
  data <- data[injury == 1]

  data[, Any_Injury_Censored := max(Censored), by = ID]

  data[, injury := NULL]
  data[, obs_id := NULL]
  data[, first_censored := NULL]

  if (specific == "delete_censored_obs") {

    data_return <- data[Censored != 1]

    data[, Any_Injury_Censored := 0]

    setcolorder(data_return, c("Injury_Length", "ID", "Censored", "Actual", "Any_Injury_Censored"))

    setkeyv(data_return, "ID")

    return(data_return)

  } else if (specific == "keep_censored_obs") {

    data_return <- data
    setcolorder(data_return, c("Injury_Length", "ID", "Censored", "Actual", "Any_Injury_Censored"))

    setkeyv(data_return, "ID")

    return(data_return)

  } else if (specific == "keep_only_max_obs") {

    data[, max_injury_length := max(Injury_Length), by = ID]

    data_return <- data[max_injury_length == Injury_Length]
    data_return[, max_injury_length := NULL]

    data_return[, delta := 1 - (1 - Censored) * Any_Injury_Censored]

    setcolorder(data_return, c("Injury_Length", "ID", "Censored", "Actual", "Any_Injury_Censored"))
    setkeyv(data_return, "ID")

    return(data_return)

  } else if (specific == "max_excess") {

    data[, max_injury_length := max(Injury_Length), by = ID]

    data_return_max <- data[max_injury_length == Injury_Length]
    data_return_max[, max_injury_length := NULL]

    data_return_max[, delta := 1 - (1 - Censored) * Any_Injury_Censored]

    data_return <- excess(data = data_return_max, ne = ne)

    data_return[, "Prop_maxima_above_thresh" := nrow(data_return) / n]


    setcolorder(data_return, c("Injury_Length", "ID", "Censored", "Actual", "Any_Injury_Censored" ))
    setkeyv(data_return, "ID")

    return(data_return)

  } else if (specific == "excess") {

    data_return <- excess(data = data, ne = ne)

    data_return[, "Prop_obs_above_thresh" := nrow(data_return) / nrow(data)]

    setcolorder(data_return, c("Injury_Length", "ID", "Censored", "Actual", "Any_Injury_Censored"))
    setkeyv(data_return, "ID")

    return(data_return)
  }
}

#' This function generates data according to \code{specific} method and then applies the maximum likelihood estimator (mle)
#' based on \code{method}.
#'
#' @param censor Integer. Where to censor observations.
#'
#' @param xi Numeric. The xi parameter in the generalized Pareto distribution.
#'
#' @param n Integer. The number of individuals.
#'
#' @param num_inj Integer. The number of injuries per individual.
#'
#' @param rate_exp The rate in the exponential distribution.
#'
#' @param ne Integer. Number of observations above the threshold. Only used when \code{specific} = 'max_excess' or
#' 'excess'.
#'
#' @param specific String. Corresponding to one of "delete_censored_obs", "keep_censored_obs", "keep_only_max_obs",
#' "max_excess" or "excess".
#'
#' @param seed Numeric. Seed to set.
#'
#' @param method String. Corresponding to one of ""MLE", "CensMLE", "MLE_delete_censored", "N", "MLE_full" or "CensMLE_full".
#'
#' @return Numeric. Corresponding to the maximum likelihood estimator based on \code{method}.
#'
#' @export
#'
#' @author Annabelle Redelmeier
#'
#' @examples
#'
#' censor <- 3; xi <- 0.5; n <- 10; num_inj <- 3; rate_exp <- 1
#' specific <- "delete_censored_obs"
#' seed <- 1
#' method <- "MLE"
#'
#' simulation(censor = censor, xi = xi, n = n, num_inj = num_inj, rate_exp = rate_exp, specific = "delete_censored_obs",
#' seed = seed, method = method)[]


simulation <- function(censor = NULL,
                       xi,
                       n,
                       num_inj,
                       rate_exp,
                       ne = NULL,
                       specific,
                       seed,
                       method){

  full <- gen_data(censor = censor,
                   xi = xi,
                   n = n,
                   num_inj = num_inj,
                   rate_exp = rate_exp,
                   ne = ne,
                   specific = specific,
                   seed = seed)

  mle(full, method = method)
}


