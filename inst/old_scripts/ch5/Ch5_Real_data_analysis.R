
library(QRM) # for peaks-of-threshold model
library(extRemes) # for exponential
library(data.table)
library(ggplot2)
library(xtable)


# For privacy reasons, no data is actually read in this script!

CdS_temp <- "read_data_here"

CdS_temp <- data.table(CdS_temp)

CdS_temp[, start := as.Data(start)]
CdS_temp[, end := as.Data(end)]
CdS_temp[, Censored := ifelse(end == "2017-01-01", 1, 0)]


# Change ids so they span the number of unique ids - use counter instead of id

CdS_UniqueId <- CdS_temp[, max := max(total_missed_perfo), by = MediCirqueID]
setkeyv(CdS_UniqueId, "MediCirqueID")
CdS_UniqueId[, counter := .I] # is .I the same as row_number()?

CdS_temp2 <- merge(CdS_UniqueId, CdS_temp, by = "MediCirqueId")
CdS_temp2 <- CdS_temp2[, c("MediCirqueId, PHCInjuryId, counter, start, end, missed_days, total_missed_perfo, Censored")]
CdS_temp2[, Injury_Length := total_missed_perfo]

# Calculate max length of injury for each injury (since sometimes they come back...?)
CdS_UniqueInjury <- CdS_temp2[, max_missed_perfo := max(total_missed_perfo), by = PHCInjuryId]
CdS_UniqueInjury[, n = .N, by = PHCInjuryID]


# this give you doubles!
CdS_temp3 <- joint(CdS_UniqueInjury, CdS_temp2,
                   by = c("PHCInjuryId" = "PHCInjuryId", "max_missed_perfo" = "Injury_Length"),
                   all.x = TRUE)

CdS <- unique(CdS_temp3, by = "PHCInjuryId")
CdS[, c("counter", "MediCirqueId", "PHCInjuryId", "start", "end", "total_missed_perfo", "max_missed_perfo", "Censored")]
CdS[, Injury_Length := max_missed_perfo]
CdS[, max_missed_pero := NULL]
setkeyv(CdS, PHCInjuryId)
CdS[, Injury_Index := .N]

# make sure the above gives the same as the original: -->
# CdS <- CdS_temp3 %>%
#   dplyr::distinct(PHCInjuryId, .keep_all = TRUE) %>%
#   dplyr::select(counter, MediCirqueId, PHCInjuryId, start, end, total_missed_perfo, max_missed_perfo, Censored) %>%
#   dplyr::mutate(Injury_Length = max_missed_perfo) %>%
#   dplyr::select(-max_missed_perfo) %>% arrange(PHCInjuryId) %>%
#   dplyr::mutate(Injury_Index = row_number())

CdS_number <- CdS[Censored == 1]

pdf("Return_PeriodAll.pdf", width = 5, height = 5)
mar.default <- c(5, 4, 4, 2) + 0.1
ggplot(data = CdS, aes(y = Injury_Length, x = Injury_Index))  + geom_jitter(col = "black") +
  xlab("Injury Id") +
  ylab("Number of Performances Missed Due to Injury") +
  geom_vline(xintercept = seq(from = 285, to = 2800, by = 285), linetype = "dotted") +
  scale_x_continuous(limits = c(-1, 2700)) +
  geom_hline(yintercept = 500) +
  theme_bw() +
  geom_text(data = data.frame(x = 0, y = 500), aes(x, y), label = "r = 500", vjust = - 0.8, hjust = -0.1)
dev.off()

# Cut off greater than 923 (3 years?)
CdS_927 <- CdS[start == "2014-01-01"]

ggplot(data = CdS_927, aes(y = Injury_Index, x = Injury_Length)) +
  geom_jitter(col = "black") +
  ylab("Injury id (number of injuries observed)") +
  xlab("Number of Performances Missed due to an Injury")  +
  geom_vline(xintercept = 926, linetype = "dotted") +
  geom_hline(yintercept = seq(from = 0, to = 20000, by = 5000)) +
  theme_bw() +
  ggtitle("5000-Observation Return Level")


pdf("Return_LevelAll.pdf", width = 5, height = 5)
mar.default <- c(5, 4, 4, 2) + 0.1
ggplot(data = CdS, aes(y = Injury_Index, x = Injury_Length)) +
  geom_jitter(col = "black") +
  ylab("Injury id (number of injuries observed)") +
  xlab("Number of Performances Missed due to an Injury")  +
  geom_vline(xintercept = 589, linetype = "dotted") +
  geom_hline(yintercept = seq(from = 0, to = 2500, by = 500)) +
  theme_bw() +
  ggtitle("500-Observation Return Level") +
  scale_y_continuous(limits = c(0, 2500))
dev.off()


###############################################################
################ Data Exploration - All Data  #################
###############################################################

p1 <- ggplot(data = CdS, aes(y = total_missed_perfo, x = counter)) +
  geom_jitter(col = "black") +
  xlab("Artist Id") +
  ylab("Number of Performances Missed Due to Injury") +
  theme_bw()

p2 <- ggplot(data = CdS, aes(total_missed_perfo)) +
  geom_histogram(binwidth = 10, col = "black", fill = "white", alpha = .2) +
  xlab("Number of Performances Missed due to Injury") +
  ylab("Count")


pdf("All.pdf", width = 5, height = 5)
mar.default <- c(5, 4, 4, 2) + 0.1
par(mar = mar.default + c(0, 5, 0, 0))
p1
dev.off()

pdf("All2.pdf", width = 5, height = 5)
mar.default <- c(5, 4, 4, 2) + 0.1
par(mar = mar.default + c(0, 5, 0, 0))
multiplot(p1, p2, cols = 2)
dev.off()

###############################################################
####### Study 1: Fit all cirque data  #########################
###############################################################

data <- CdS

pdf("ME.pdf", width = 5, height = 5)
mar.default <- c(5, 4, 4, 2) + 0.1
MEplot(data = data$Injury_Length, main = "", xaxt = "n")
axis(1, at = c(110, 330, 460, 600, 800), label = TRUE)
abline(v = 110, col = "grey")
abline(v = 330, col = "grey")
abline(v = 460, col = "grey")
dev.off()

pdf("xi.pdf", width = 5, height = 5)
mar.default <- c(5, 4, 4, 2) + 0.1
xiplot(data = data$Injury_Length, models = 30., start = 100, end = 2000., reverse = TRUE,
       ci = 0.95)
abline(h = 0.045, col = "grey", lty = 2)
abline(v = - 940, col = "grey")
abline(v = - 550, col = "grey")
dev.off()

threshold <- 100
data <- CdS
data_u <- data[Injury_Length > 100]
data_u[, Injury_Length := Injury_Length - 100]

data_u_cens <- data_u[Censored == 1]

CdS_ties <- data_u[, n := .N, by = "Injury_Length"]

method <- c("MLE_full","CensMLE_full")

met_all <- matrix(NA, nrow = 2, ncol = 4)
colnames(met_all) <- c("MLE", "se", "CensMLE", "se")
rownames(met_all) <- c("xi", "beta")


for (i in 1:length(method)) {
  met_all[1, 2 * i - 1] <- mle(data = data, threshold = threshold, method = method[i])
  met_all[1, i + i] <-  mle(data = data, threshold = threshold, method = method[i])$par.ses[1]
  met_all[2, 2 * i - 1] <-  mle(data = data, threshold = threshold, method = method[i],
                             information = "observed")$par.ests[2]
  met_all[2, i + i] <-  mle(data = data, threshold = threshold, method = method[i],
                         information = "observed")$par.ses[2]
}
met_all


###############################################################
################ Calculate Quantiles  #########################
###############################################################

out <- fit.GPD(data = data$total_missed_perfo, threshold = 100)

RiskMeasures(out, c(0.95, 0.99))

# compare with exponential

gpd_extRemes <- extRemes::fevd(x = data$total_missed_perfo, threshold = 100, type = "GP")
gpd_extRemes_exp <- extRemes::fevd(x = data$total_missed_perfo, threshold = 100,
                                   type = "Exponential")

extRemes::lr.test(gpd_extRemes_exp, gpd_extRemes)

###############################################################
##### Under Assumption that Injury Excess Exponential  ########
###############################################################

threshold <- 100

data <- CdS
data_u <- data %>%
  dplyr::filter(Injury_Length > threshold) %>%
  dplyr::mutate(Injury_Length = Injury_Length - threshold)

beta_reg <- as.numeric(extRemes::fevd(x = data$total_missed_perfo, threshold = 100,
                                      type = "Exponential")$results$par)

r_1 <- 200
r_2 <- 500
r_3 <- 1000

u <- 100
Nu <- nrow(data_u)
n <- nrow(data)
N <- 2356

#P(X > xk) for exponential

Prob_pareto_hat1 <- 1 / ((Nu / n) * exp(- (r_1 - u) / beta_reg))
Prob_pareto_hat2 <- 1 / ((Nu / n) * exp(- (r_2 - u) / beta_reg))
Prob_pareto_hat3 <- 1 / ((Nu / n) * exp(- (r_3 - u) / beta_reg))


RP <- matrix(NA, nrow = 1, ncol = 3)
colnames(RP) <- c(r_1, r_2, r_3)
rownames(RP) <- c("1/P(X>r)")
RP[1, ] <- c(Prob_pareto_hat1, Prob_pareto_hat2, Prob_pareto_hat3)

xtable(RP, digits = c(1, 3, 3, 3))

# the k-level return level is
k_1 <- 50
k_2 <- 500
k_3 <- 5000


zk_pareto_hat1 <- beta_reg * log(k_1 * Nu / n) + u
zk_pareto_hat2 <- beta_reg * log(k_2 * Nu / n) + u
zk_pareto_hat3 <- beta_reg * log(k_3 * Nu / n) + u


RL <- matrix(NA, nrow = 1, ncol = 3)
colnames(RL) <- c(k_1, k_2, k_3)
rownames(RL) <- c("r_k")
RL[1, ] <- c(zk_pareto_hat1, zk_pareto_hat2, zk_pareto_hat3)

xtable(RL, digits = c(1, 3, 3, 3))

# Plot return level
CdS_950 <- CdS[start <= "2014-01-01"]

ggplot(data = CdS_950, aes(y = Injury_Index, x = Injury_Length)) +
  geom_jitter(col = "black") +
  ylab("Injury id (number of injuries observed)") +
  xlab("Number of Performances Missed Due to an Injury") +
  geom_vline(xintercept = 950, linetype = "dotted") +
  geom_hline(yintercept = c(5000, 10000)) + theme_bw() +
  ggtitle("5000-Observation Return Level")

pdf("Return_LevelAll.pdf", width = 5, height = 5)
mar.default <- c(5, 4, 4, 2) + 0.1
ggplot(data = CdS_950, aes(y = Injury_Index, x = Injury_Length)) +
  geom_jitter(col = "black") +
  ylab("Injury Id") +
  xlab("Number of Performances Missed Due to Injury")  +
  geom_vline(xintercept = 950, linetype = "dotted") +
  geom_hline(yintercept = c(5000, 10000)) +
  theme_bw() +
  geom_text(data = data.frame(x = 1000, y = 5000), aes(x, y), label = "k=5000", hjust = - 0.5, vjust = -1)
dev.off()


# QQ Plot (I changed the qqline function a little to work with the exponential distribution)

qqline <- function (y, datax = FALSE, probs = c(0.25, 0.75), qtype = 7, ...) {
  stopifnot(length(probs) == 2, is.function(distribution))
  y <- quantile(y, probs, names = FALSE, type = qtype, na.rm = TRUE)
  x <- qexp(probs, rate = 1 / beta_reg)
  if (datax) {
    slope <- diff(x) / diff(y)
    int <- x[1L] - slope * y[1L]
  }
  else {
    slope <- diff(y)/diff(x)
    int <- y[1L] - slope * x[1L]
  }
  abline(int, slope, ...)
}

data <- data_u$Injury_Length

pdf("QQ_All.pdf", width = 5, height = 5)
mar.default <- c(5, 4, 4, 2) + 0.1
qqplot(x = qexp(p = ppoints(100),
              rate = 1 / beta_reg), y = data,
       main = "Exponential Q-Q Plot",
       xlab = "Theoretical Quantiles",
       ylab = "Sample Quantiles")
qqline(data)
dev.off()

###############################################################
################ Data Exploration - Just Maxima  ##############
###############################################################

ggplot(data = CdS, aes(y = total_missed_perfo, x = counter)) +
  geom_jitter(col = "black") +
  ggtitle("Total Missed Performance per Artist") +
  xlab("Artist Id") +
  ylab("Total Missed Performance") +
  coord_cartesian(xlim = c(50, 150), ylim = c(0, 750))


p19 <- ggplot(data = CdS, aes(y = total_missed_perfo, x = counter))  +
  geom_jitter(col = "black") +
  ggtitle("Total Missed Performance per Artist") +
  xlab("Artist Id") +
  ylab("Total Missed Performance") +
  coord_cartesian(xlim = c(19, 25), ylim = c(0, 750)) +
  geom_vline(xintercept = 19:25, linetype = "dotted")

p175 <- ggplot(data = CdS, aes(y = total_missed_perfo, x = counter))  +
  geom_jitter(col = "black") +
  ggtitle("Total Missed Performance per Artist") +
  xlab("Artist Id") +
  ylab("Total Missed Performance") +
  coord_cartesian(xlim = c(175, 200), ylim = c(0, 750))


###############################################################
################### Just Maxima  ##############################
###############################################################

cirque_max <- CdS[, max_missed_perfo := max(total_missed_perfo), by = "counter"]
cirque_max[, min_start := min(start), by = "counter"]

CdS.2 <- merge(CdS, cirque_max, by = "counter", by.x = TRUE)
CdS.2[, Longest_Injury := ifelse(max_missed_perfo == total_missed_perfo, 1, 0)]
CdS.2[, Longest_Injury := as.factor(Longest_Injury)]

p4 <- ggplot(data = cirque_max, aes(max_missed_perfo)) +
  geom_histogram(binwidth = 10, col = "black", fill = "white", alpha = .2) +
  coord_cartesian(xlim = c(0, 800)) +
  ggtitle("Maximum Missed Performance per Artist") +
  xlab("Maximum Number Missed Performance") +
  ylab("Count")


###############################################################
####### Study 2: Fit only maxima  #############################
###############################################################

tmp <- CdS[, Injury_Length_max := max(total_missed_perfo), by = "MediCirqueID"]

data_max_double <- merge(tmp, CdS,
                         by = c("MediCirqueId" = "MediCirqueId", "Injury_Length_max" = "total_missed_perfo"),
                         all.x = TRUE)

data_max <- unique(data_max_double, by = "MediCirqueId")
# ? same as below ?
# data_max <- dplyr::distinct(data_max_double, MediCirqueId, .keep_all = TRUE)

pdf("Return_PeriodMax.pdf", width = 5, height = 5)
mar.default <- c(5, 4, 4, 2) + 0.1
ggplot(data = data_max, aes(y = Injury_Length, x = counter))  +
  geom_jitter(col = "black") +
  xlab("Artist id") +
  ylab("Number of Performances Missed due to Injury")  +
  geom_vline(xintercept = seq(from = 0, to = 2441, by = 706), linetype = "dotted") +
  geom_hline(yintercept = 1000) +
  theme_bw()
dev.off()

pdf("Return_PeriodMax.pdf", width = 5, height = 5)
mar.default <- c(5, 4, 4, 2) + 0.1
ggplot(data = data_max, aes(y = Injury_Length, x = counter))  +
  geom_jitter(col = "black") +
  xlab("Artist id") +
  ylab("Number of Performances Missed due to Injury")  +
  geom_vline(xintercept = seq(from = 0, to = 240, by = 40), linetype = "dotted") +
  scale_x_continuous(limits = c(0, 240)) +
  geom_hline(yintercept = 500)
dev.off()

pdf("ME_max.pdf", width = 5, height = 5)
mar.default <- c(5, 4, 4, 2) + 0.1
MEplot(data = data_max$Injury_Length_max, main = "", xaxt = "n")
axis(1, at = c(110, 250, 370, 600, 800), label = TRUE)
#abline(v = 80, col = "grey")
abline(v = 110, col = "grey")
abline(v = 370, col = "grey")
abline(v = 430, col = "grey")
dev.off()

pdf("xi_max.pdf", width = 5, height = 5)
mar.default <- c(5, 4, 4, 2) + 0.1
xiplot(data = data_max$Injury_Length_max, models = 30., start = 15., end = 1500., reverse = TRUE,
       ci = 0.95)
abline(h = 0, col = "grey")
abline(v = -731, col = "grey")
abline(v = -450, col = "grey")
dev.off()

threshold <- 100

data_max_u <- data_max[Injury_Length_max > threshold]
data_max_u[, Injury_Length := Injury_Length - threshold]

data <- data_max

method <- c("MLE_full", "CensMLE_full")

met_all_max <- matrix(NA, nrow = 2, ncol = 4)
colnames(met_all_max) <- c("MLE", "se", "CensMLE", "se")
rownames(met_all_max) <- c("xi", "beta")

# this is my function mle that's why you need to input all the data

for (i in 1:length(method)) {
  met_all_max[1, 2 * i - 1] <- mle(data = data, threshold = threshold, method = method[i],
                                 information = "observed")$par.ests[1]
  met_all_max[1, i + i] <- mle(data = data, threshold = threshold, method = method[i],
                             information = "observed")$par.ses[1]
  met_all_max[2, 2 * i - 1] <- mle(data = data, threshold = threshold, method = method[i],
                                 information = "observed")$par.ests[2]
  met_all_max[2, i + i] <- mle(data = data, threshold = threshold, method = method[i],
                             information = "observed")$par.ses[2]
}
met_all_max

# compare with exponential

gpd_extRemes_max <- extRemes::fevd(x = data$Injury_Length, threshold = threshold, type = "GP")
gpd_extRemes_exp_max <- extRemes::fevd(x = data$Injury_Length, threshold = threshold, type = "Exponential")

lr.test(gpd_extRemes_exp_max, gpd_extRemes_max)

###############################################################
########## Return level under exponential model  ##############
###############################################################

beta_reg <- as.numeric(extRemes::fevd(x = data$Injury_Length, threshold = 100, type = "Exponential")$results$par)

r_1 <- 200
r_2 <- 500
r_3 <- 1000

u <- 100
Nu <- nrow(data_max_u)
n <- nrow(data_max)

# Return Periods

#P(X > xk) for exponential
Prob_pareto_hat1 <- (Nu / n) * exp(- (r_1 - u) / beta_reg)
Prob_pareto_hat2 <- (Nu / n) * exp(- (r_2 - u) / beta_reg)
Prob_pareto_hat3 <- (Nu / n) * exp(- (r_3 - u) / beta_reg)

RP <- matrix(NA, nrow = 1, ncol = 3)
colnames(RP) <- c(r_1, r_2, r_3)
rownames(RP) <- c("1/P(X>r)")
RP[1, ] <- c(1 / Prob_pareto_hat1, 1 / Prob_pareto_hat2, 1 / Prob_pareto_hat3)

xtable(RP, digits = c(3, 3, 3, 3))


# Return Levels
k_1 <- 50
k_2 <- 500
k_3 <- 1500

# the k-level return level is
zk_pareto_hat1 <- beta_reg * log(k_1 * Nu / n) + u
zk_pareto_hat2 <- beta_reg * log(k_2 * Nu / n) + u
zk_pareto_hat3 <- beta_reg * log(k_3 * Nu / n) + u

RL <- matrix(NA, nrow = 1, ncol = 3)
colnames(RL) <- c(k_1, k_2, k_3)
rownames(RL) <- c("r_k")
RL[1, ] <- c(zk_pareto_hat1, zk_pareto_hat2, zk_pareto_hat3)

xtable(RL, digits = c(3, 3, 3, 3))


# QQ Plots

# This is data after subtracting the threshold 100
data <- data_max_u$Injury_Length

pdf("QQ_Max.pdf", width = 5, height = 5)
mar.default <- c(5, 4, 4, 2) + 0.1
qqplot(x = qexp(p = ppoints(100), rate = 1 / beta_reg), y = data, main = "Exponential Q-Q Plot",
       xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(data)
dev.off()

###############################################################
########## Pick threshold and decluster after  ################
###############################################################


# 95 percent quantile
# quantile(CdS$total_missed_perfo, .90)
# quantile(CdS$total_missed_perfo, .95) = 88
# quantile(CdS$total_missed_perfo, .99) = 337

data88 <- data[total_missed_perfo > 88]
data88[, n := .N, by = "MediCirqueId"][n > 1]
# ? same ?
# dplyr::filter(total_missed_perfo > 88) %>%
# group_by(MediCirqueId) %>%
# dplyr::summarise(n = n()) %>%
# dplyr::filter(n > 1)

# 208 artists have more than 2 injuries greater than 81
# 208 / 958 = 22%

data337 <- data[total_missed_perfo > 337]
data337[, n := .N, by = "MediCirqueId"][n > 1]
# 36 artists have more than 2 injuries greater than 323
# 36/194 = 19%


CdS.2 <- merge(CdS, cirque_max, by = "counter", all.x = TRUE)
CdS.2[, Longest_Injury := ifelse(max_missed_perfo == total_missed_perfo, 1, 0)]
CdS.2[, Longest_Injury := as.factor(Longest_Injury)]

pmax <- ggplot(data = CdS.2, aes(y = total_missed_perfo, x = counter, color = Longest_Injury)) +
  geom_jitter() +
  scale_y_continuous(breaks = sort(c(0, 81, 323, 500, 1000))) +
  xlab("Artist Id") +
  ylab("Total Missed Performances")


pmax_u <- ggplot(data = CdS.2, aes(y = total_missed_perfo, x = counter, color = Longest_Injury)) +
  geom_jitter() +
  scale_y_continuous(breaks = sort(c(0, 88, 337, 500, 1000))) +
  xlab("Artist Id") +
  ylab("Number of Performances Missed Due to Injury") +
  geom_hline(yintercept = c(88, 337), col = "black") +
  theme_bw()

pdf("pmax_u.pdf", width = 7, height = 4)
mar.default <- c(5, 4, 4, 2) + 0.1
pmax_u
dev.off()

ggplot(data = CdS.2, aes(y = total_missed_perfo, x = counter, color = Longest_Injury)) +
  geom_jitter() +
  geom_hline(yintercept = c(81, 323), col = "black") +
  scale_y_continuous(breaks = sort(c(0, 81, 323, 500, 1000))) +
  xlab("Artist Id") +
  ylab("Total Missed Performances") +
  coord_cartesian(xlim = c(108, 120)) +
  geom_vline(xintercept = 0:125, linetype = "dotted")


CdS.81 <- CdS[total_missed_perfo > 81]
CdS.81[, n := .N, by = "counter"][n >= 2]

CdS.12 <- CdS[total_missed_perfo == 12]

ggplot(data = CdS.12, aes(y = total_missed_perfo, x = start)) +
  geom_point() +
  geom_hline(yintercept = c(81), col = "black") +
  scale_y_continuous(breaks = sort(c(0, 81, 100, 200, 300))) +
  geom_vline(xintercept = c(as.Date("2008-06-8"),
                            as.Date("2008-09-15"),
                            as.Date("2011-05-01"),
                            as.Date("2012-10-15"),
                            as.Date("2013-03-31"),
                            as.Date("2014-03-02")),
             linetype = "dotted") +
  xlab("Date")

CdS.323 <- CdS[total_missed_perfo > 323]
CdS.323[, n := .N, by = "counter"][n >= 2]



###############################################################
#################### Other Plots  #############################
###############################################################

# Can we do it as a function of when the injury occured?

cirque_when <- merge(cirque_max, CdS, by = c("MediCirqueId" = "MediCirqueId",
                                                 "max_missed_perfo" = "total_missed_perfo"),
                     all.x = TRUE)

cirque_when[, diff := start - min_start]
cirque_when[, c("MediCirqueId", "max_missed_perfo", "min_start", "start", "diff")]

p5 <- ggplot(data = cirque_when, aes(y = max_missed_perfo, diff / 365)) +
  geom_jitter(col = "blue") +
  xlab("Time Since First Injury Date (Years)") +
  ylab("Maximum Number Missed Performances") +
  ggtitle("Time Since First Injury") +
  theme(text = element_text(size = 22))

# Number of injuries per artist

cirque_num <- CdS[, n = .N, by = "MediCirqueId"]

p6 <- ggplot(data = cirque_num, aes(numb)) +
  geom_histogram(binwidth = 1, col = "blue", fill = "light blue", alpha = 0.2) +
  xlab("Number of Injuries") +
  ggtitle("Number of Completely-Out Injuries per Artist") +
  ylab("Count")

# ggplot(data = cirque_all, aes(y = total_missed_perfo, x = start, group = MediCirqueId)) + geom_jitter()

CdS[, month := format(start, "%Y/%m")]

cirque_month <- CdS[, mean := mean(total_missed_perfo), by = "month"]


cirque_year <- merge(cirque_max, CdS, by = c("MediCirqueId" = "MediCirqueId",
                                                 "max_missed_perfo" = "total_missed_perfo"),
                     all.x = TRUE)

cirque_year[, c("MediCirqueId", "max_missed_perfo", "start")]
cirque_year[, year := format(start, "%Y")]
cirque_year[, year := as.numeric(year)]


cirque_year_num <- cirque_year[, n := N, by = "year"]


pdf("cirque_year.pdf", width = 5, height = 5)
mar.default <- c(5, 4, 4, 2) + 0.1
barplot(cirque_year_num$n, xlab = "Year",
        names.arg = c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016"))
abline(h = 2356)
mtext(text = "Number of Injuries", side = 2, line = 3, las = 0)
dev.off()


p7 <- ggplot(data = cirque_year, aes(y = total_missed_perfo, x = year)) +
  geom_point(color = "black") +
  geom_smooth(method = "lm", color = "blue", se = FALSE) +
  ggtitle("Max Missed Performance per Year") +
  xlab("Year") +
  ylab("Maximum Number Missed Performances")
