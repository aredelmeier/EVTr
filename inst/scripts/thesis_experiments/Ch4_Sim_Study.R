
library(ggplot2)
library(data.table)
library(EVTr)

## To make Figure 4-3: Simulation 3

n.sim <- 10
n <- 100
num_inj <- 10
xi <- c(0.1, 0.3, 0.5)
rate_exp <- 1
ne <- 5
method <- c("MLE")

results_MLE <- matrix(NA, nrow = n.sim, ncol = length(xi))

results <- list()
for (met in method) {
  for (i in 1:n.sim) {
    row_MLE <- NULL
    for (x in xi) {
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

df[, xi_num := as.numeric(as.character(xi))]
ggplot(df, aes(x = method, y = value - xi_num, fill = xi)) +
  geom_boxplot() +
  xlab("") +
  ylab("Observed Bias") +
  geom_hline(yintercept = 0, col = "red")


## To make Table 4-1: Simulation 3

df0 <- df[, c("value", "ID", "xi_num")]


df_agg <- aggregate(. ~ ID, df0, function(x) c(bias = mean(x), mean = mean(x), sd = sd(x)))
df_agg <- data.table(df_agg)
df_agg[, value.bias := value.bias - xi_num.bias]

df_agg[, xi_num.bias := NULL]
df_agg[, xi_num.mean := NULL]
df_agg[, xi_num.sd := NULL]


df1 <- df[, c("ID", "method", "xi")]
df1 <- unique(df1)

df_all <- merge(df_agg, df1, by = "ID")


# Simulation 4 - Figure 4-4
n.sim <- 100 # should be 1000
n <- 100 # should be 1000
num_inj <- 60
xi <- 0.5
censor <- c(36, 60, 90)
rate_exp <- 1
ne <- 5 # should be 50?
method <- c()


results_Norm_list <- list()

results_Norm <- matrix(NA, nrow = n.sim, ncol = length(censor))

results <- list()

for (i in 1:n.sim) {

  Norm_Value <- NULL
  for (cens in censor) {
    full_Norm_value <- gen_data(cens,
                                xi,
                                n,
                                num_inj,
                                rate_exp,
                                ne,
                                specific = "max_excess",
                                seed = i)
    Norm_Value <- c(Norm_Value, full_Norm_value$Prop_maxima_censored[1] * 100)
  }
  results_Norm[i, ] <- Norm_Value

}
results_Norm_list[["Norm_Value_max"]] <- data.table(results_Norm)

results_MLE_dt <- rbindlist(results_Norm_list, idcol = TRUE)

results_MLE_dt <- cbind(results_MLE_dt, data_manu = "max_excess")
colnames(results_MLE_dt) <- c("method", as.character(censor), "data_manu")

results_MLE_dt[, n := 1:.N, by = c("method", "data_manu")]

df <- data.table::melt(results_MLE_dt, id.vars = c("method", "n", "data_manu"), variable.name = "censor")
setkeyv(df, "method")

# Figure 4-4
ggplot(df, aes(x = method, y = value, fill = censor)) +
  geom_boxplot() +
  xlab("End-of-Followup Time") +
  ylab("Percentage of Censored Observation")




# Simulation 4 - Table 4-2

n.sim <- 100 # should be 1000
n <- 100 # should be 1000
num_inj <- 60
xi <- 0.5
censor <- c(36, 60, 90)
rate_exp <- 1
ne <- 5 # should be 50?
method <- c("MLE", "CensMLE")

results_Norm_list <- list()
results_Cens_list <- list()

results_Norm <- matrix(NA, nrow = n.sim, ncol = length(censor))
results_Cens <- matrix(NA, nrow = n.sim, ncol = length(censor))

results <- list()

for (i in 1:n.sim) {
  Norm_Value <- NULL
  Cens_Value <- NULL
  for (cens in censor) {
    full_Norm_value <- gen_data(cens,
                                xi,
                                n,
                                num_inj,
                                rate_exp,
                                ne,
                                specific = "max_excess",
                                seed = i)


    Norm_Value <- c(Norm_Value, mle(full_Norm_value, method = "MLE"))

    #
    full_CensValue <- gen_data(cens,
                               xi,
                               n,
                               num_inj,
                               rate_exp,
                               ne,
                               specific = "max_excess",
                               seed = i)




    Cens_Value <- c(Cens_Value, mle(full_CensValue, method = "CensMLE"))


  }
  results_Norm[i, ] <- Norm_Value
  results_Cens[i, ] <- Cens_Value

}
results_Norm_list[["MLE"]] <- data.table(results_Norm)
results_Cens_list[["CensMLE"]] <- data.table(results_Cens)

results_Norm_dt <- rbindlist(results_Norm_list, idcol = TRUE)
results_Norm_dt <- cbind(results_Norm_dt, data_manu = "max_excess")

results_Cens_dt <- rbindlist(results_Cens_list, idcol = TRUE)
results_Cens_dt <- cbind(results_Cens_dt, data_manu = "max_excess")

results_MLE_dt <- rbind(results_Norm_dt, results_Cens_dt)

colnames(results_MLE_dt) <- c("method", as.character(censor), "data_manu")
results_MLE_dt[, n := 1:.N, by = c("method", "data_manu")]

df <- data.table::melt(results_MLE_dt, id.vars = c("method", "n", "data_manu"), variable.name = "censor")
setkeyv(df, "method")

ggplot(df, aes(x = censor, y = value, fill = method)) +
  geom_boxplot() +
  xlab("End-of_Followup Time") +
  ylab("Average Injury Length") +
  geom_hline(yintercept = 0, col = "red")



# Table 4-2

df[, ID := .GRP, by = .(method, data_manu, censor)]

df0 <- df[, c("value", "ID")]

df_agg <- aggregate(. ~ ID, df0, function(x) c(mean = mean(x), bias = mean(x) - 0.5, sd = sd(x)))

df1 <- df[, c("ID", "method", "data_manu", "censor")]
df1 <- unique(df1)

df_all <- merge(df_agg, df1, by = "ID")
df_all <- data.table(df_all)

tmp_non_cens_like <- df_all[method == "MLE"][, c("censor", "value.bias", "value.sd", "value.mean")]

tmp_cens_like <- df_all[method == "CensMLE"][, c("censor", "value.bias", "value.sd", "value.mean")]

cbind(tmp_non_cens_like, tmp_cens_like[, censor := NULL])

