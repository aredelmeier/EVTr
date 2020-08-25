library(EVTr)
library(QRM)
library(data.table)

###############################################################
####### Study 1: without censoring  ###########################
###############################################################

n.sim = 10
n = 100
num_inj = 30
xi = c(0.1, 0.5, 1.5)
rate_exp = 0.1
ne <- 5 # this is used in "excess"
method = c("MLE", "CensMLE")

results_MLE <- matrix(NA, nrow = n.sim, ncol = length(xi))

results <- list()
for(met in method) {
  for(i in 1:n.sim) {
    row_MLE = NULL
    for(x in xi) {
      row_MLE <- c(row_MLE, simulation(censor = NULL,
                                       xi = x,
                                       n = n,
                                       num_inj = num_inj,
                                       rate_exp = rate_exp,
                                       ne = ne,
                                       method = met,
                                       specific = "max_excess",
                                       seed = i))


    }
    results_MLE[i, ] <- row_MLE

  }
  results[[met]] <- data.table(results_MLE)

}

results_MLE_dt <- rbindlist(results, idcol = TRUE)

colnames(results_MLE_dt) <- c("method", as.character(xi))
results_MLE_dt[, n := 1:.N, by = method]

df <- data.table::melt(results_MLE_dt, id.vars = c("method", "n"), variable.name = "xi")
setkeyv(df, "method")


df[, ID := .GRP, by = .(method, xi)]

df0 <- df[, c("value", "ID")]

df_agg <- aggregate(. ~ ID, df0, function(x) c(mean = mean(x), sd = sd(x)))

df1 <- df[, c("ID", "method", "xi")]
df1 <- unique(df1)

df_all <- merge(df_agg, df1, by = "ID")


pdf("Sim_March_No_Cens.pdf", width = 5, height = 5)
mar.default <- c(5, 4, 4, 2) + 0.1
par(mar = mar.default + c(0, 5, 0, 0))
boxplot(I(value - as.numeric(as.character(df$xi))) ~ xi, data = df)
abline(h = 0, col = "red")
mtext(text = "True xi", side = 1, line = 4, las = 1)
mtext(text = "Observed Bias", side = 2, line = 4, las = 0)
dev.off()


###############################################################
####### Study 1.5: without censoring (April)  #################
###############################################################


n.sim = 10
n = 100
num_inj = 10
xi = c(0.1, 0.5, 1.5)
rate_exp = 1
ne <- 5 # this is used in "excess"
method = c("MLE", "CensMLE")

results_MLE <- matrix(NA, nrow = n.sim, ncol = length(xi))

results <- list()
for(met in method) {
  for(i in 1:n.sim) {
    row_MLE = NULL
    for(x in xi) {
      row_MLE <- c(row_MLE, simulation(censor = NULL,
                                       xi = x,
                                       n = n,
                                       num_inj = num_inj,
                                       rate_exp = rate_exp,
                                       ne = ne,
                                       method = met,
                                       specific = "max_excess",
                                       seed = i))


    }
    results_MLE[i, ] <- row_MLE

  }
  results[[met]] <- data.table(results_MLE)

}

results_MLE_dt <- rbindlist(results, idcol = TRUE)

colnames(results_MLE_dt) <- c("method", as.character(xi))
results_MLE_dt[, n := 1:.N, by = method]

df <- data.table::melt(results_MLE_dt, id.vars = c("method", "n"), variable.name = "xi")
setkeyv(df, "method")


df[, ID := .GRP, by = .(method, xi)]

df0 <- df[, c("value", "ID")]

df_agg <- aggregate(. ~ ID, df0, function(x) c(mean = mean(x), sd = sd(x)))

df1 <- df[, c("ID", "method", "xi")]
df1 <- unique(df1)

df_all <- merge(df_agg, df1, by = "ID")


pdf("Sim_March_No_Cens.pdf", width = 5, height = 5)
mar.default <- c(5, 4, 4, 2) + 0.1
par(mar = mar.default + c(0, 5, 0, 0))
boxplot(I(value - as.numeric(as.character(df$xi))) ~ xi, data = df)
abline(h = 0, col = "red")
mtext(text = "True xi", side = 1, line = 4, las = 1)
mtext(text = "Observed Bias", side = 2, line = 4, las = 0)
dev.off()


###############################################################
#### Study 1.5: without censoring (April - change xi)  ########
###############################################################

n.sim = 10
n = 100
num_inj = 10
xi = c(0.1, 0.3, 0.5)
rate_exp = 1
ne <- 5 # this is used in "excess"
method = c("MLE", "CensMLE")

results_MLE <- matrix(NA, nrow = n.sim, ncol = length(xi))

results <- list()
for(met in method) {
  for(i in 1:n.sim) {
    row_MLE = NULL
    for(x in xi) {
      row_MLE <- c(row_MLE, simulation(censor = NULL,
                                       xi = x,
                                       n = n,
                                       num_inj = num_inj,
                                       rate_exp = rate_exp,
                                       ne = ne,
                                       method = met,
                                       specific = "max_excess",
                                       seed = i))


    }
    results_MLE[i, ] <- row_MLE

  }
  results[[met]] <- data.table(results_MLE)

}



results_MLE_dt <- rbindlist(results, idcol = TRUE)

colnames(results_MLE_dt) <- c("method", as.character(xi))
results_MLE_dt[, n := 1:.N, by = method]

df <- data.table::melt(results_MLE_dt, id.vars = c("method", "n"), variable.name = "xi")
setkeyv(df, "method")

df[, ID := .GRP, by = .(method, xi)]

df0 <- df[, c("value", "ID")]

df_agg <- aggregate(. ~ ID, df0, function(x) c(mean = mean(x), sd = sd(x)))

df1 <- df[, c("ID", "method", "xi")]
df1 <- unique(df1)

df_all <- merge(df_agg, df1, by = "ID")


pdf("Sim_March_No_Cens.pdf", width = 5, height = 5)
mar.default <- c(5, 4, 4, 2) + 0.1
par(mar = mar.default + c(0, 5, 0, 0))
boxplot(I(value - as.numeric(as.character(df$xi))) ~ xi, data = df)
abline(h = 0, col = "red")
mtext(text = "True xi", side = 1, line = 4, las = 1)
mtext(text = "Observed Bias", side = 2, line = 4, las = 0)
dev.off()


###############################################################
####### Study 2: with censoring  ##############################
###############################################################


n.sim = 10
n = 10
num_inj = 60
xi = 0.5
censor <- c(36, 60, 90)
rate_exp = 1
ne <- 5 # this is used in "excess"
method = c("MLE", "CensMLE")

results_MLE <- matrix(NA, nrow = n.sim, ncol = length(censor))

results <- list()
for(met in method) {
  for(i in 1:n.sim) {
    row_MLE = NULL
    for(cens in censor) {


      row_MLE <- c(row_MLE, simulation(censor = cens,
                                       xi = xi,
                                       n = n,
                                       num_inj = num_inj,
                                       rate_exp = rate_exp,
                                       ne = ne,
                                       method = met,
                                       specific = "max_excess",
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

df[, ID := .GRP, by = .(method, censor)]

df0 <- df[, c("value", "ID")]

df_agg <- aggregate(. ~ ID, df0, function(x) c(mean = mean(x), sd = sd(x)))

df1 <- df[, c("ID", "method", "censor")]
df1 <- unique(df1)

df_all <- merge(df_agg, df1, by = "ID")


pdf("Sim_March_No_Cens.pdf", width = 5, height = 5)
mar.default <- c(5, 4, 4, 2) + 0.1
par(mar = mar.default + c(0, 5, 0, 0))
boxplot(I(value - as.numeric(as.character(df$censor))) ~ censor, data = df)
abline(h = 0, col = "red")
mtext(text = "True censor", side = 1, line = 4, las = 1)
mtext(text="Observed Bias", side = 2, line = 4, las = 0)
dev.off()





###############################################################
####### Study 3: count censored  ##############################
###############################################################


n.sim = 10
n = 100
num_inj = 60
xi = 0.5
censor <- c(36, 60, 90)
rate_exp = 1
ne <- 5 # this is used in "excess"
method = c("N")

results_MLE <- matrix(NA, nrow = n.sim, ncol = length(censor))

results <- list()
for(met in method) {
  for(i in 1:n.sim) {
    row_MLE = NULL
    for(cens in censor) {


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


results_MLE_dt <- rbindlist(results, idcol = TRUE)

colnames(results_MLE_dt) <- c("method", as.character(censor))
results_MLE_dt[, n := 1:.N, by = method]

df <- data.table::melt(results_MLE_dt, id.vars = c("method", "n"), variable.name = "censor")
setkeyv(df, "method")

pdf("Sim_March_N.pdf", width = 5, height = 5)
mar.default <- c(5, 4, 4, 2) + 0.1
par(mar = mar.default + c(0, 5, 0, 0))
boxplot(I(value) ~ censor, data = df, xaxt = "n", yaxt = "n")
axis(side=1, at= 1:3, labels= c("36", "60", "120"))
axis(side=2, at= c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7), labels= c("20", "30", "40", "50", "60", "70" ))
mtext(text="End-of-Followup Time", side = 1, line = 4, las = 1)
mtext(text="Percentage of Censored Observations", side = 2, line = 4, las = 0)
dev.off()




