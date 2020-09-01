
# This R code calculates the mean injury length for non-censored injuries vs censored injuries

library(EVTr)
library(QRM)
library(data.table)


################################################################
####################### All Observations #######################
################################################################


################################################################
######################### Theta = 1 ############################
################################################################

n.sim <- 10
n <- 10
num_inj <- 30
xi <- 1
censor <- c(3, 6, 12, 36)
rate_exp <- 1
method <- c("Norm_Value", "Cens_Value")


results_Norm_list <- list()
results_Cens_list <- list()

results_Norm <- matrix(NA, nrow = n.sim, ncol = length(censor))
results_Cens <- matrix(NA, nrow = n.sim, ncol = length(censor))

results <- list()

for (i in 1:n.sim) {
  for (cens in censor) {
    Norm_Value <- NULL
    Cens_Value <- NULL

    # Norm_value --> data_nc <- data[data[, "Censored"] == 0]
    full_Norm_value <- gen_data(cens,
                                xi,
                                n,
                                num_inj,
                                rate_exp,
                                ne,
                                specific = "delete_censored_obs",
                                seed = i)


    Norm_Value <- c(Norm_Value, mean(full_Norm_value$Injury_Length))

    # Cens_Value --> # data_c <- data[data[, "Censored"] == 1]
    # mean(data$Actual)
    full_CensValue <- gen_data(cens,
                               xi,
                               n,
                               num_inj,
                               rate_exp,
                               ne,
                               specific = "keep_censored_obs",
                               seed = i)


    full_CensValue0 <- full_CensValue[Censored == 1]

    Cens_Value <- c(Cens_Value, mean(full_CensValue0[["Actual"]]))


  }
  results_Norm[i, ] <- Norm_Value
  results_Cens[i, ] <- Cens_Value

}
results_Norm_list[["Norm_Value"]] <- data.table(results_Norm)
results_Cens_list[["Cens_Value"]] <- data.table(results_Cens)


results_Norm_dt <- rbindlist(results_Norm_list, idcol = TRUE)
results_Norm_dt <- cbind(results_Norm_dt, data_manu = "delete_censored_obs")

results_Cens_dt <- rbindlist(results_Cens_list, idcol = TRUE)
results_Cens_dt <- cbind(results_Cens_dt, data_manu = "keep_censored_obs")

results_MLE_dt <- rbind(results_Norm_dt, results_Cens_dt)

colnames(results_MLE_dt) <- c("method", as.character(censor), "data_manu")
results_MLE_dt[, n := 1:.N, by = c("method", "data_manu")]

df <- data.table::melt(results_MLE_dt, id.vars = c("method", "n", "data_manu"), variable.name = "censor")
setkeyv(df, "method")


################################################################
######################### Theta = 1 April ######################
################################################################


n.sim <- 10
n <- 10
num_inj <- 60
xi <- 1
censor <- c(36, 60, 90)
rate_exp <- 1
method <- c("Norm_Value", "Cens_Value")

results_Norm_list <- list()
results_Cens_list <- list()

results_Norm <- matrix(NA, nrow = n.sim, ncol = length(censor))
results_Cens <- matrix(NA, nrow = n.sim, ncol = length(censor))

results <- list()

for (i in 1:n.sim) {
  for (cens in censor) {
    Norm_Value <- NULL
    Cens_Value <- NULL

    # Norm_value --> data_nc <- data[data[, "Censored"] == 0]
    full_Norm_value <- gen_data(cens,
                                xi,
                                n,
                                num_inj,
                                rate_exp,
                                ne,
                                specific = "delete_censored_obs",
                                seed = i)


    Norm_Value <- c(Norm_Value, mean(full_Norm_value$Injury_Length))

    # Cens_Value --> # data_c <- data[data[, "Censored"] == 1]
    # mean(data$Actual)
    full_CensValue <- gen_data(cens,
                               xi,
                               n,
                               num_inj,
                               rate_exp,
                               ne,
                               specific = "keep_censored_obs",
                               seed = i)


    full_CensValue0 <- full_CensValue[Censored == 1]

    Cens_Value <- c(Cens_Value, mean(full_CensValue0[["Actual"]]))


  }
  results_Norm[i, ] <- Norm_Value
  results_Cens[i, ] <- Cens_Value

}
results_Norm_list[["Norm_Value"]] <- data.table(results_Norm)
results_Cens_list[["Cens_Value"]] <- data.table(results_Cens)


results_Norm_dt <- rbindlist(results_Norm_list, idcol = TRUE)
results_Norm_dt <- cbind(results_Norm_dt, data_manu = "delete_censored_obs")

results_Cens_dt <- rbindlist(results_Cens_list, idcol = TRUE)
results_Cens_dt <- cbind(results_Cens_dt, data_manu = "keep_censored_obs")

results_MLE_dt <- rbind(results_Norm_dt, results_Cens_dt)

colnames(results_MLE_dt) <- c("method", as.character(censor), "data_manu")
results_MLE_dt[, n := 1:.N, by = c("method", "data_manu")]

df <- data.table::melt(results_MLE_dt, id.vars = c("method", "n", "data_manu"), variable.name = "censor")
setkeyv(df, "method")


df_Normal36 <- df[method == "Norm_Value"][censor == "36"]
df_Normal36 <- df[method == "Cens_Value"][censor == "36"]

pdf("Sim3_hist_all.pdf", width = 5, height = 5)
mar.default <- c(5, 4, 4, 2) + 0.1
par(mar = mar.default + c(0, 5, 0, 0))
hist(df_Normal36$value, main = "theta = 1, non-censored obs, c = 36", xlab = "")
mtext(text = "Average Injury Length", side = 1, line = 4, las = 1)
dev.off()

pdf("Sim3_hist_all2.pdf", width = 5, height = 5)
mar.default <- c(5, 4, 4, 2) + 0.1
par(mar = mar.default + c(0, 5, 0, 0))
hist(df_Censor36$value, breaks = seq(0, (max(df_Censor36$value) + 50), by = 40),
     main = "theta = 1, censored obs, c = 36", xlab = "")
mtext(text = "Average Injury Length", side = 1, line = 4, las = 1)
dev.off()

pdf("Sim3_1.pdf", width = 6, height = 5)
mar.default <- c(5, 4, 4, 2) + 0.1
par(mar = mar.default + c(0, 5, 0, 0))
boxplot(I(value) ~ method + censor, data = df, xaxt = "n")
abline(v = 2.5, col = "grey")
abline(v = 4.5, col = "grey")
axis(side = 1, at =  1:6, label = c("36", "36", "60", "60", "90", "90"))
mtext(text = "End-of-Followup Time", side = 1, line = 4, las = 1)
mtext(text = "Average Injury Length", side = 2, line = 4, las = 0)
dev.off()

pdf("Sim3_1_2.pdf", width = 5, height = 5)
mar.default <- c(5, 4, 4, 2) + 0.1
par(mar = mar.default + c(0, 5, 0, 0))
boxplot(I(value) ~ method + censor, data = df, xaxt = "n", main = "theta = 1", ylim = c(0, 50))
abline(v = 2.5, col = "grey")
abline(v = 4.5, col = "grey")
axis(side = 1, at = 1:6, label = c("12", "12", "24", "24", "36", "36"))
mtext(text = "Censoring Time", side = 1, line = 4, las = 1)
mtext(text = "Average Injury Length", side = 2, line = 4, las = 0)
dev.off()


################################################################
######################### Theta = 0.2 ##########################
################################################################


n.sim <- 10
n <- 10
num_inj <- 30
xi <- 0.5
censor <- c(36, 90, 120)
rate_exp <- 0.2
method <- c("Norm_Value", "Cens_Value")

results_Norm_list <- list()
results_Cens_list <- list()

results_Norm <- matrix(NA, nrow = n.sim, ncol = length(censor))
results_Cens <- matrix(NA, nrow = n.sim, ncol = length(censor))

results <- list()

for (i in 1:n.sim) {
  for (cens in censor) {
    Norm_Value <- NULL
    Cens_Value <- NULL

    # Norm_value --> data_nc <- data[data[, "Censored"] == 0]
    full_Norm_value <- gen_data(cens,
                                xi,
                                n,
                                num_inj,
                                rate_exp,
                                ne,
                                specific = "delete_censored_obs",
                                seed = i)


    Norm_Value <- c(Norm_Value, mean(full_Norm_value$Injury_Length))

    # Cens_Value --> # data_c <- data[data[, "Censored"] == 1]
    # mean(data$Actual)
    full_CensValue <- gen_data(cens,
                               xi,
                               n,
                               num_inj,
                               rate_exp,
                               ne,
                               specific = "keep_censored_obs",
                               seed = i)


    full_CensValue0 <- full_CensValue[Censored == 1]

    Cens_Value <- c(Cens_Value, mean(full_CensValue0[["Actual"]]))


  }
  results_Norm[i, ] <- Norm_Value
  results_Cens[i, ] <- Cens_Value

}
results_Norm_list[["Norm_Value"]] <- data.table(results_Norm)
results_Cens_list[["Cens_Value"]] <- data.table(results_Cens)


results_Norm_dt <- rbindlist(results_Norm_list, idcol = TRUE)
results_Norm_dt <- cbind(results_Norm_dt, data_manu = "delete_censored_obs")

results_Cens_dt <- rbindlist(results_Cens_list, idcol = TRUE)
results_Cens_dt <- cbind(results_Cens_dt, data_manu = "keep_censored_obs")

results_MLE_dt <- rbind(results_Norm_dt, results_Cens_dt)

colnames(results_MLE_dt) <- c("method", as.character(censor), "data_manu")
results_MLE_dt[, n := 1:.N, by = c("method", "data_manu")]

df <- data.table::melt(results_MLE_dt, id.vars = c("method", "n", "data_manu"), variable.name = "censor")
setkeyv(df, "method")


df_Normal120 <- df[method == "Norm_Value"][censor == "120"]
df_Normal120 <- df[method == "Cens_Value"][censor == "120"]

pdf("Sim3_hist_all_2.pdf", width = 5, height = 5)
mar.default <- c(5, 4, 4, 2) + 0.1
par(mar = mar.default + c(0, 5, 0, 0))
hist(df2_Normal120$value, main = "theta = 0.2, non-censored obs, c = 120", xlab = "")
mtext(text = "Average Injury Length", side = 1, line = 4, las = 1)
dev.off()

pdf("Sim3_hist_all_2_2.pdf", width = 5, height = 5)
mar.default <- c(5, 4, 4, 2) + 0.1
par(mar = mar.default + c(0, 5, 0, 0))
hist(df2_Censor120$value, breaks = seq(0, 350, by = 2), xlim = c(0, 250),
     main = "theta = 0.2, censored obs, c = 120", xlab = "")
mtext(text = "Average Injury Length", side = 1, line = 4, las = 1)
dev.off()

pdf("Sim3_2.pdf", width = 5, height = 5)
mar.default <- c(5, 4, 4, 2) + 0.1
par(mar = mar.default + c(0, 5, 0, 0))
boxplot(I(value) ~ method + censor, data = df2,  xaxt = "n", main = "theta = 0.2")
abline(v = 2.5, col = "grey")
abline(v = 4.5, col = "grey")
axis(side = 1, at = 1:6, label = c("36", "36", "90", "90", "120", "120"))
mtext(text = "Censoring Time", side = 1, line = 4, las = 1)
mtext(text = "Average Injury Length", side = 2, line = 4, las = 0)
dev.off()

pdf("Sim3_2_2.pdf", width = 5, height = 5)
mar.default <- c(5, 4, 4, 2) + 0.1
par(mar = mar.default + c(0, 5, 0, 0))
boxplot(I(value)~method + censor, data = df2,  xaxt = "n", main = "theta = 0.2", ylim = c(0, 50))
abline(v = 2.5, col = "grey")
abline(v = 4.5, col = "grey")
axis(side = 1, at = 1:6, label =  c("36", "36", "90", "90", "120", "120"))
mtext(text = "Censoring Time", side = 1, line = 4, las = 1)
mtext(text = "Average Injury Length", side = 2, line = 4, las = 0)
dev.off()

################################################################
####################### Just Maxima ############################
################################################################


################################################################
######################### Theta = 1 ############################
################################################################

n.sim <- 10
n <- 10
num_inj <- 30
xi <- 0.5
censor <- c(24, 36)
rate_exp <- 1
ne <- 50
method <- c("Norm_Value_max", "Cens_Value")


results_Norm_list <- list()
results_Cens_list <- list()

results_Norm <- matrix(NA, nrow = n.sim, ncol = length(censor))
results_Cens <- matrix(NA, nrow = n.sim, ncol = length(censor))

results <- list()

for (i in 1:n.sim) {
  for (cens in censor) {
    Norm_Value <- NULL
    Cens_Value <- NULL

    # Norm_Value_max --> data_nc <- data[data[, "Censored"] == 0]
    # mean(data_nc$Injury_Length_before)
    full_Norm_value <- gen_data(cens,
                                xi,
                                n,
                                num_inj,
                                rate_exp,
                                ne,
                                specific = "excess",
                                seed = i)


    Norm_Value <- c(Norm_Value, mean(full_Norm_value$Injury_Length_before))

    # Cens_Value --> # data_c <- data[data[, "Censored"] == 1]
    # mean(data$Actual)
    full_CensValue <- gen_data(cens,
                               xi,
                               n,
                               num_inj,
                               rate_exp,
                               ne,
                               specific = "keep_censored_obs",
                               seed = i)


    full_CensValue0 <- full_CensValue[Censored == 1]

    Cens_Value <- c(Cens_Value, mean(full_CensValue0[["Actual"]]))


  }
  results_Norm[i, ] <- Norm_Value
  results_Cens[i, ] <- Cens_Value

}
results_Norm_list[["Norm_Value_max"]] <- data.table(results_Norm)
results_Cens_list[["Cens_Value"]] <- data.table(results_Cens)


results_Norm_dt <- rbindlist(results_Norm_list, idcol = TRUE)
results_Norm_dt <- cbind(results_Norm_dt, data_manu = "delete_censored_obs")

results_Cens_dt <- rbindlist(results_Cens_list, idcol = TRUE)
results_Cens_dt <- cbind(results_Cens_dt, data_manu = "keep_censored_obs")

results_MLE_dt <- rbind(results_Norm_dt, results_Cens_dt)

colnames(results_MLE_dt) <- c("method", as.character(censor), "data_manu")
results_MLE_dt[, n := 1:.N, by = c("method", "data_manu")]

df <- data.table::melt(results_MLE_dt, id.vars = c("method", "n", "data_manu"), variable.name = "censor")
setkeyv(df, "method")


df_Normal36 <- df[method == "Norm_Value_max"][censor == "36"]
df_Censor36 <- df[method == "Cens_Value"][censor == "36"]


pdf("Sim3_hist1.pdf", width = 5, height = 5)
mar.default <- c(5, 4, 4, 2) + 0.1
par(mar = mar.default + c(0, 5, 0, 0))
hist(df_Normal36$value, main = "theta = 1, non-censored maxima,
     c = 36", xlab = "")
mtext(text = "Average Injury Length", side = 1, line = 4, las = 1)
dev.off()

pdf("Sim3_hist2.pdf", width = 5, height = 5)
mar.default <- c(5, 4, 4, 2) + 0.1
par(mar = mar.default + c(0, 5, 0, 0))
hist(df_Censor36$value, breaks = seq(0, 4000, by = 20), main = "theta = 1, censored maxima,
     c = 36", xlab = "")
mtext(text = "Average Injury Length", side = 1, line = 4, las = 1)
dev.off()

pdf("Sim3_max_1.pdf", width = 5, height = 5)
mar.default <- c(5, 4, 4, 2) + 0.1
par(mar = mar.default + c(0, 5, 0, 0))
boxplot(I(value - 0.5) ~ method + censor, data = df, xaxt = "n", main = "theta = 1 *")
abline(v = 2.5, col = "grey")
axis(side = 1, at = 1:4, label = c("24", "24", "36", "36"))
mtext(text = "Censoring Time", side = 1, line = 4, las = 1)
mtext(text = "Average Injury Length", side = 2, line = 4, las = 0)
dev.off()

pdf("Sim3_max_1_2.pdf", width = 5, height = 5)
mar.default <- c(5, 4, 4, 2) + 0.1
par(mar = mar.default + c(0, 5, 0, 0))
boxplot(I(value - 0.5) ~ method + censor, data = df, xaxt = "n", main = "theta = 1 *", ylim = c(0, 400))
abline(v = 2.5, col = "grey")
axis(side = 1, at = 1:4, label = c("24", "24", "36", "36"))
mtext(text = "Censoring Time", side = 1, line = 4, las = 1)
mtext(text = "Average Injury Length", side = 2, line = 4, las = 0)
dev.off()

################################################################
######################### Theta = 0.2 ##########################
################################################################


n.sim <- 10
n <- 10
num_inj <- 30
xi <- 0.5
censor <- c(90, 120)
rate_exp <- 0.2
ne <- 50
method <- c("Norm_Value_max", "Cens_Value")

results_Norm_list <- list()
results_Cens_list <- list()

results_Norm <- matrix(NA, nrow = n.sim, ncol = length(censor))
results_Cens <- matrix(NA, nrow = n.sim, ncol = length(censor))

results <- list()

for (i in 1:n.sim) {
  for (cens in censor) {
    Norm_Value <- NULL
    Cens_Value <- NULL

    # Norm_Value_max --> data_nc <- data[data[, "Censored"] == 0]
    # mean(data_nc$Injury_Length_before)
    full_Norm_value <- gen_data(cens,
                                xi,
                                n,
                                num_inj,
                                rate_exp,
                                ne,
                                specific = "excess",
                                seed = i)


    Norm_Value <- c(Norm_Value, mean(full_Norm_value$Injury_Length_before))

    # Cens_Value --> # data_c <- data[data[, "Censored"] == 1]
    # mean(data$Actual)
    full_CensValue <- gen_data(cens,
                               xi,
                               n,
                               num_inj,
                               rate_exp,
                               ne,
                               specific = "keep_censored_obs",
                               seed = i)


    full_CensValue0 <- full_CensValue[Censored == 1]

    Cens_Value <- c(Cens_Value, mean(full_CensValue0[["Actual"]]))


  }
  results_Norm[i, ] <- Norm_Value
  results_Cens[i, ] <- Cens_Value

}
results_Norm_list[["Norm_Value_max"]] <- data.table(results_Norm)
results_Cens_list[["Cens_Value"]] <- data.table(results_Cens)


results_Norm_dt <- rbindlist(results_Norm_list, idcol = TRUE)
results_Norm_dt <- cbind(results_Norm_dt, data_manu = "delete_censored_obs")

results_Cens_dt <- rbindlist(results_Cens_list, idcol = TRUE)
results_Cens_dt <- cbind(results_Cens_dt, data_manu = "keep_censored_obs")

results_MLE_dt <- rbind(results_Norm_dt, results_Cens_dt)

colnames(results_MLE_dt) <- c("method", as.character(censor), "data_manu")
results_MLE_dt[, n := 1:.N, by = c("method", "data_manu")]

df <- data.table::melt(results_MLE_dt, id.vars = c("method", "n", "data_manu"), variable.name = "censor")
setkeyv(df, "method")

df_Normal120 <- df[method == "Norm_Value_max"][censor == "120"]
df_Censor120 <- df[method == "Cens_Value"][censor == "120"]

pdf("Sim3_hist_2.pdf", width = 5, height = 5)
mar.default <- c(5, 4, 4, 2) + 0.1
par(mar = mar.default + c(0, 5, 0, 0))
hist(df_Normal120$value, main = "theta = 0.2, non-censored maxima,
     c = 120", xlab = "")
mtext(text = "Average Injury Length", side = 1, line = 4, las = 1)
dev.off()

pdf("Sim3_hist_2_2.pdf", width = 5, height = 5)
mar.default <- c(5, 4, 4, 2) + 0.1
par(mar = mar.default + c(0, 5, 0, 0))
hist(df_Censor120$value, breaks = seq(0, 6000, by = 50), main = "theta = 0.2, censored maxima,
     c = 120", xlab = "")
mtext(text = "Average Injury Length", side = 1, line = 4, las = 1)
dev.off()

pdf("Sim3_max_2.pdf", width = 5, height = 5)
mar.default <- c(5, 4, 4, 2) + 0.1
par(mar = mar.default + c(0, 5, 0, 0))
boxplot(I(value - 0.5) ~ method + censor, data = df, main = "theta = 0.2 *", xaxt = "n")
abline(v = 2.5, col = "grey")
axis(side = 1, at = 1:4, label = c("90", "90", "120", "120"))
mtext(text = "Censoring Time", side = 1, line = 4, las = 1)
mtext(text = "Average Injury Length", side = 2, line = 4, las = 0)
dev.off()

pdf("Sim3_max_2_2.pdf", width = 5, height = 5)
mar.default <- c(5, 4, 4, 2) + 0.1
par(mar = mar.default + c(0, 5, 0, 0))
boxplot(I(value - 0.5) ~ method + censor, data = df, main = "theta = 0.2", xaxt = "n", ylim = c(0, 400))
abline(v = 2.5, col = "grey")
axis(side = 1, at = 1:4, label = c("90", "90", "120", "120"))
mtext(text = "Censoring Time", side = 1, line = 4, las = 1)
mtext(text = "Average Injury Length", side = 2, line = 4, las = 0)
dev.off()
