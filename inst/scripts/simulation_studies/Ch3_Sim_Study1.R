library(EVTr)
library(QRM)
library(data.table)

###############################################################
######### Study 1: Delete Censored Obs  #######################
###############################################################

n.sim <- 10
n <- 10
num_inj <- 30
xi <- 0.5
censor <- c(3, 6, 12, 36)
rate_exp <- 1
method <- c("MLE", "CensMLE")

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

df <- data.table::melt(results_MLE_dt, id.vars = c("method", "n"), variable.name = "censor")
setkeyv(df, "method")

pdf("Sim2_delete_March.pdf", width = 5, height = 5)
mar.default <- c(5, 4, 4, 2) + 0.1
par(mar = mar.default + c(0, 5, 0, 0))
boxplot(I(value - 0.5) ~ censor + method, data = df, xaxt = "n", main = "MLE/CensMLE")
abline(h = 0, col = "red")
abline(v = 4.5, col = "grey")
axis(side = 1, at = 1:8, label = c("3", "6", "12", "36", "3", "6", "12", "36"))
mtext(text = "Censoring Time", side = 1, line = 4, las = 1)
mtext(text = "Observed Bias", side = 2, line = 4, las = 0)
dev.off()


###############################################################
######### Study 2: Keep Censored Obs  #########################
###############################################################

n.sim <- 10
n <- 10
num_inj <- 30
xi <- 0.5
censor <- c(3, 6, 12, 36)
rate_exp <- 1
method <- c("MLE", "CensMLE")
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

df <- data.table::melt(results_MLE_dt, id.vars = c("method", "n"), variable.name = "censor")
setkeyv(df, "method")


pdf("Sim2_keep_March.pdf", width = 5, height = 5)
mar.default <- c(5, 4, 4, 2) + 0.1
par(mar = mar.default + c(0, 5, 0, 0))
boxplot(I(value - 0.5) ~ censor + method, data = df, xaxt = "n")
abline(h = 0, col = "red")
abline(v = 4.5, col = "grey")
axis(side = 1, at = 1:8, labels = c("3", "6", "12", "36", "3", "6", "12", "36"))
mtext(text = "Censoring Time", side = 1, line = 4, las = 0)
mtext(text = "Observed Bias", side = 2, line = 4, las = 0)
dev.off()

###############################################################
####### Study 2 Extra: Both Keep and Delete (all data)  #######
###############################################################

n.sim <- 10
n <- 1000
num_inj <- 30
xi <- 0.5
censor <- c(3, 6, 12, 36)
rate_exp <- 1
method <- c("MLE", "CensMLE")


results_delete_censored_list <- list()
results_keep_censored_list <- list()

results_delete_censored_obs <- matrix(NA, nrow = n.sim, ncol = length(censor))
results_keep_censored_obs <- matrix(NA, nrow = n.sim, ncol = length(censor))

for (met in method) {
  row_MLE <- NULL

  for (i in 1:n.sim) {
    for (cens in censor) {
      delete_censored_obs <- c(row_MLE, simulation(censor = cens,
                                       xi = xi,
                                       n = n,
                                       num_inj = num_inj,
                                       rate_exp = rate_exp,
                                       method = met,
                                       specific = "delete_censored_obs",
                                       seed = i))

      keep_censored_obs <- c(row_MLE, simulation(censor = cens,
                                                   xi = xi,
                                                   n = n,
                                                   num_inj = num_inj,
                                                   rate_exp = rate_exp,
                                                   method = met,
                                                   specific = "keep_censored_obs",
                                                   seed = i))
    }
    results_delete_censored_obs[i, ] <- delete_censored_obs
    results_keep_censored_obs[i, ] <- keep_censored_obs
  }



  results_delete_censored_list[[met]] <- data.table(results_delete_censored_obs)
  results_keep_censored_list[[met]] <- data.table(results_keep_censored_obs)

}


results_MLE_delete <- rbindlist(results_delete_censored_list, idcol = TRUE)
results_MLE_delete <- cbind(results_MLE_delete, data_manu = "delete_censored_obs")

results_MLE_keep <- rbindlist(results_keep_censored_list, idcol = TRUE)
results_MLE_keep <- cbind(results_MLE_keep, data_manu = "keep_censored_obs")

results_MLE_dt <- rbind(results_MLE_delete, results_MLE_keep)

colnames(results_MLE_dt) <- c("method", as.character(censor), "data_manu")
results_MLE_dt[, n := 1:.N, by = c("method", "data_manu")]

df <- data.table::melt(results_MLE_dt, id.vars = c("method", "n", "data_manu"), variable.name = "censor")
setkeyv(df, "method")

pdf("Sim2_all.pdf", width = 7, height = 5)
mar.default <- c(5, 4, 4, 2) + 0.1
par(mar = mar.default + c(0, 5, 0, 0))
boxplot(I(value - 0.5) ~ censor + method + data_manu, data = df, xaxt = "n",
        main = paste("Delete MLE/", "Keep MLE/", "Delete CensMLE/", "Keep CensMLE"))
abline(h = 0, col = "red")
abline(v = 4.5, col = "grey")
abline(v = 8.5, col = "grey")
abline(v = 12.5, col = "grey")
axis(side = 1, at = 1:(length(censor) * 4), labels = c(rep(censor, 4)))
mtext(text = "End-Of-Followup Time", side = 1, line = 4, las = 0)
mtext(text = "Observed Bias", side = 2, line = 4, las = 0)
dev.off()


df[, ID := .GRP, by = .(method, data_manu, censor)]

df0 <- df[, c("value", "ID")]

df_agg <- aggregate(. ~ ID, df0, function(x) c(mean = mean(x), sd = sd(x)))

df1 <- df[, c("ID", "method", "data_manu", "censor")]
df1 <- unique(df1)

df_all <- merge(df_agg, df1, by = "ID")


###############################################################
####### Study 2.5 Extra: Both Keep and Delete (all data)  #######
###############################################################

n.sim <- 10
n <- 10
num_inj <- 60
xi <- 0.5
censor <- c(30, 60, 90)
rate_exp <- 1
method <- c("MLE", "CensMLE")

results_delete_censored_list <- list()
results_keep_censored_list <- list()

results_delete_censored_obs <- matrix(NA, nrow = n.sim, ncol = length(censor))
results_keep_censored_obs <- matrix(NA, nrow = n.sim, ncol = length(censor))

for (met in method) {
  row_MLE <- NULL

  for (i in 1:n.sim) {
    for (cens in censor) {
      delete_censored_obs <- c(row_MLE, simulation(censor = cens,
                                                   xi = xi,
                                                   n = n,
                                                   num_inj = num_inj,
                                                   rate_exp = rate_exp,
                                                   method = met,
                                                   specific = "delete_censored_obs",
                                                   seed = i))

      keep_censored_obs <- c(row_MLE, simulation(censor = cens,
                                                 xi = xi,
                                                 n = n,
                                                 num_inj = num_inj,
                                                 rate_exp = rate_exp,
                                                 method = met,
                                                 specific = "keep_censored_obs",
                                                 seed = i))
    }
    results_delete_censored_obs[i, ] <- delete_censored_obs
    results_keep_censored_obs[i, ] <- keep_censored_obs
  }



  results_delete_censored_list[[met]] <- data.table(results_delete_censored_obs)
  results_keep_censored_list[[met]] <- data.table(results_keep_censored_obs)

}


results_MLE_delete <- rbindlist(results_delete_censored_list, idcol = TRUE)
results_MLE_delete <- cbind(results_MLE_delete, data_manu = "delete_censored_obs")

results_MLE_keep <- rbindlist(results_keep_censored_list, idcol = TRUE)
results_MLE_keep <- cbind(results_MLE_keep, data_manu = "keep_censored_obs")

results_MLE_dt <- rbind(results_MLE_delete, results_MLE_keep)

colnames(results_MLE_dt) <- c("method", as.character(censor), "data_manu")
results_MLE_dt[, n := 1:.N, by = c("method", "data_manu")]

df <- data.table::melt(results_MLE_dt, id.vars = c("method", "n", "data_manu"), variable.name = "censor")
setkeyv(df, "method")

pdf("Sim2_all.pdf", width = 7, height = 5)
mar.default <- c(5, 4, 4, 2) + 0.1
par(mar = mar.default + c(0, 5, 0, 0))
boxplot(I(value - 0.5) ~ censor + method + data_manu, data = df, xaxt = "n",
        main = paste("Delete MLE", "Keep MLE", "Delete CensMLE", "Keep CensMLE"))
abline(h = 0, col = "red")
abline(v = 4.5, col = "grey")
abline(v = 8.5, col = "grey")
abline(v = 12.5, col = "grey")
axis(side = 1, at = 1:(length(censor) * 4), labels = c(rep(censor, 4)))
mtext(text = "End-Of-Followup Time", side = 1, line = 4, las = 0)
mtext(text = "Observed Bias", side = 2, line = 4, las = 0)
dev.off()


df[, ID := .GRP, by = .(method, data_manu, censor)]

df0 <- df[, c("value", "ID")]

df_agg <- aggregate(. ~ ID, df0, function(x) c(mean = mean(x), sd = sd(x)))

df1 <- df[, c("ID", "method", "data_manu", "censor")]
df1 <- unique(df1)

df_all <- merge(df_agg, df1, by = "ID")

###############################################################
####### Study 3: Delete Censored Obs (Only max)  ##############
###############################################################

n.sim <- 10
n <- 100
num_inj <- 30
xi <- 0.5
censor <- c(3, 6, 12, 36)
rate_exp <- 1
method <- c("MLE", "CensMLE")

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
                                       specific = "keep_only_max_obs",
                                       seed = i))
    }
    results_MLE[i, ] <- row_MLE

  }
  results[[met]] <- data.table(results_MLE)
}


results_MLE_dt <- rbindlist(results, idcol = TRUE)

colnames(results_MLE_dt) <- c("method", as.character(censor))
results_MLE_dt[, n := 1:.N, by = method]

df <- data.table::melt(results_MLE_dt, id.vars = c("method", "n"), variable.name = "censor")
setkeyv(df, "method")

pdf("Sim2_delete_March.pdf", width = 5, height = 5)
mar.default <- c(5, 4, 4, 2) + 0.1
par(mar = mar.default + c(0, 5, 0, 0))
boxplot(I(value - 0.5) ~ censor + method, data = df, xaxt = "n", main = "MLE/CensMLE")
abline(h = 0, col = "red")
abline(v = 4.5, col = "grey")
axis(side = 1, at = 1:8, label = c("3", "6", "12", "36", "3", "6", "12", "36"))
mtext(text = "Censoring Time", side = 1, line = 4, las = 1)
mtext(text = "Observed Bias", side = 2, line = 4, las = 0)
dev.off()

###############################################################
####### Study 4: Keep Censored Obs (Only max)  ################
###############################################################

n.sim <- 10
n <- 100
num_inj <- 30
xi <- 0.5
censor <- c(30, 60, 90, 120)
rate_exp <- 0.2
method <- c("MLE", "CensMLE")

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
                                       specific = "keep_only_max_obs",
                                       seed = i))
    }
    results_MLE[i, ] <- row_MLE

  }
  results[[met]] <- data.table(results_MLE)
}


results_MLE_dt <- rbindlist(results, idcol = TRUE)

colnames(results_MLE_dt) <- c("method", as.character(censor))
results_MLE_dt[, n := 1:.N, by = method]

df <- data.table::melt(results_MLE_dt, id.vars = c("method", "n"), variable.name = "censor")
setkeyv(df, "method")


pdf("Cens_keep02_March.pdf", width = 5, height = 5)
mar.default <- c(5, 4, 4, 2) + 0.1
par(mar = mar.default + c(0, 1, 0, 0))
boxplot(I(value - 0.5) ~ censor + method, data = df, xaxt = "n", main = "MLE/CensMLE")
abline(h = 0, col = "red")
abline(v = 4.5, col = "grey")
axis(side = 1, at = 1:8, labels = c("36", "60", "90", "120", "36", "60", "90", "120"))
mtext(text = "Censoring Time", side = 1, line = 4, las = 1)
mtext(text = "Observed Bias", side = 2, line = 4, las = 0)
dev.off()
