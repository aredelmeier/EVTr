
library(ggplot2)
library(data.table)
library(EVTr)

# Simulation 3.4 (thesis)

n.sim <- 10 # thesis switch this to 1000
n <- 10 # thesis switch this to 1000
num_inj <- 60
xi <- 0.5
censor <- c(3, 6, 12, 36)
rate_exp <- 1


# 1. Treat censored observations as not censored

method <- c("MLE")

results_MLE <- matrix(NA, nrow = n.sim, ncol = length(censor))

results <- list()
for (met in method) {
  for (i in 1:n.sim) {
    row_MLE <- NULL
    row_CensMLE <- NULL

    for (cens in censor) {

      row_MLE <- c(row_MLE, simulation(censor = cens,
                                       xi = xi,
                                       n = n,
                                       num_inj = num_inj,
                                       rate_exp = rate_exp,
                                       method = met,
                                       specific = "keep_censored_obs",
                                       seed = i))
    }
    results_MLE[i, ] <- row_MLE

  }
  results[[met]] <- data.table(results_MLE)
}


results_MLE_dt <- rbindlist(results, idcol = TRUE)

colnames(results_MLE_dt) <- c("method", as.character(censor))
results_MLE_dt[, n := 1:.N, by = method]

df1 <- data.table::melt(results_MLE_dt, id.vars = c("method", "n"), variable.name = "censor")
setkeyv(df1, "method")


# 2. Delete censored observations

method <- c("MLE")

results_MLE <- matrix(NA, nrow = n.sim, ncol = length(censor))

results <- list()
for (met in method) {
  for (i in 1:n.sim) {
    row_MLE <- NULL
    row_CensMLE <- NULL

    for (cens in censor) {

      row_MLE <- c(row_MLE, simulation(censor = cens,
                                       xi = xi,
                                       n = n,
                                       num_inj = num_inj,
                                       rate_exp = rate_exp,
                                       method = met,
                                       specific = "delete_censored_obs",
                                       seed = i))
    }
    results_MLE[i, ] <- row_MLE

  }
  results[[met]] <- data.table(results_MLE)
}


results_MLE_dt <- rbindlist(results, idcol = TRUE)

colnames(results_MLE_dt) <- c("method", as.character(censor))
results_MLE_dt[, n := 1:.N, by = method]

df2 <- data.table::melt(results_MLE_dt, id.vars = c("method", "n"), variable.name = "censor")
setkeyv(df2, "method")

# 3.

method <- c("CensMLE")

results_MLE <- matrix(NA, nrow = n.sim, ncol = length(censor))

results <- list()
for (met in method) {
  for (i in 1:n.sim) {
    row_MLE <- NULL
    row_CensMLE <- NULL

    for (cens in censor) {

      row_MLE <- c(row_MLE, simulation(censor = cens,
                                       xi = xi,
                                       n = n,
                                       num_inj = num_inj,
                                       rate_exp = rate_exp,
                                       method = met,
                                       specific = "keep_censored_obs",
                                       seed = i))
    }
    results_MLE[i, ] <- row_MLE

  }
  results[[met]] <- data.table(results_MLE)
}


results_MLE_dt <- rbindlist(results, idcol = TRUE)

colnames(results_MLE_dt) <- c("method", as.character(censor))
results_MLE_dt[, n := 1:.N, by = method]

df3 <- data.table::melt(results_MLE_dt, id.vars = c("method", "n"), variable.name = "censor")
setkeyv(df3, "method")


# To get Table 3-1: Simulation 2

df1 <- cbind(df1, data_manu = "keep_censored_obs")
df2 <- cbind(df2, data_manu = "delete_censored_obs")
df3 <- cbind(df3, data_manu = "keep_censored_obs")
df <- rbind(df1, df2, df3)

df[, ID := .GRP, by = .(method, data_manu, censor)]

df0 <- df[, c("value", "ID")]

df_agg <- aggregate(. ~ ID, df0, function(x) c(mean = mean(x), sd = sd(x)))

df1 <- df[, c("ID", "method", "data_manu", "censor")]
df1 <- unique(df1)

df_all <- merge(df_agg, df1, by = "ID")
df_all

# some how also need to come up with bias, sd, MSE?

# To get Figure 3-4: Simulation 2

# this looks slightly different then the thesis plot - but only because
# it was easier to use ggplot2 this time around
df[, method_datamanu := paste(method, data_manu, sep = "-")]

ggplot(df, aes(x = method_datamanu, y = value - 0.5, fill = censor)) +
  geom_boxplot() +
  xlab("mle method + data manipulation") +
  ylab("bias")

# Figure 3-3

# Make a figure such that for different censoring times we have the
# average length of the censored observations
# and the average length of the non-censored observations

method <- c("N")

results_MLE <- matrix(NA, nrow = n.sim, ncol = length(censor))


results <- list()
for (met in method) {
  for (i in 1:n.sim) {
    row_MLE <- NULL
    for (cens in censor) {


      data <- gen_data(censor = cens,
                       xi = xi,
                       n = n,
                       num_inj = num_inj,
                       rate_exp = rate_exp,
                       ne = ne,
                       specific = "max_excess",
                       seed = i)

      row_MLE <- c(row_MLE, data$Percent_of_obs_above_thresh[1])

    }
    results_MLE[i, ] <- row_MLE

  }
  results[[met]] <- data.table(results_MLE)

}





