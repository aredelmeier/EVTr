
library(QRM) # for peaks-of-threshold model
library(extRemes) # for exponential

library(data.table)
library(ggplot2)
library(xtable)


# For privacy reasons, no data is actually read in this script!

# Figure 5-1
p1 <- ggplot(data = CdS, aes(y = total_missed_perfo, x = ID)) +
  geom_jitter(col = "black") +
  xlab("Artist Id") +
  ylab("Number of Performances Missed Due to Injury") +
  theme_bw()


pdf("All.pdf", width = 5, height = 5)
mar.default <- c(5, 4, 4, 2) + 0.1
par(mar = mar.default + c(0, 5, 0, 0))
p1
dev.off()

# Figure 5-2

pdf("ME.pdf", width = 5, height = 5)
mar.default <- c(5, 4, 4, 2) + 0.1
QRM::MEplot(data = CdS$Injury_Length, main = "", xaxt = "n")
axis(1, at = c(110, 330, 460, 600, 800), label = TRUE)
abline(v = 110, col = "grey")
abline(v = 330, col = "grey")
abline(v = 460, col = "grey")
dev.off()

# Figure 5-3


pdf("xi.pdf", width = 5, height = 5)
mar.default <- c(5, 4, 4, 2) + 0.1
QRM::xiplot(data = data$Injury_Length, models = 30., start = 100, end = 2000., reverse = TRUE,
       ci = 0.95)
abline(h = 0.045, col = "grey", lty = 2)
abline(v = - 940, col = "grey")
abline(v = - 550, col = "grey")
dev.off()


# Table 5-1


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


## Table 5-2
threshold <- 100

data <- CdS
data_u <- data[Injury_Length > threshold)]
data_u[, Injury_Length := Injury_Length - threshold]

beta_reg <- as.numeric(extRemes::fevd(x = data$total_missed_perfo, threshold = 100,
                                      type = "Exponential")$results$par)

r_1 <- 200
r_2 <- 500
r_3 <- 1000

u <- 100
Nu <- nrow(data_u)
n <- nrow(data)
N <- 2356

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

# Figure 5-5


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


# Figure 5-6

cirque_max <- CdS[, max_missed_perfo := max(total_missed_perfo), by = "counter"]
cirque_max[, min_start := min(start), by = "counter"]

CdS.2 <- merge(CdS, cirque_max, by = "counter", all.x = TRUE)
CdS.2[, Longest_Injury := ifelse(max_missed_perfo == total_missed_perfo, 1, 0)]
CdS.2[, Longest_Injury := as.factor(Longest_Injury)]


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

# Figure 5-7

tmp <- CdS[, Injury_Length_max := max(total_missed_perfo), by = "MediCirqueID"]

data_max_double <- merge(tmp, CdS,
                         by = c("MediCirqueId" = "MediCirqueId", "Injury_Length_max" = "total_missed_perfo"),
                         all.x = TRUE)

data_max <- unique(data_max_double, by = "MediCirqueId")

pdf("ME_max.pdf", width = 5, height = 5)
mar.default <- c(5, 4, 4, 2) + 0.1
MEplot(data = data_max$Injury_Length_max, main = "", xaxt = "n")
axis(1, at = c(110, 250, 370, 600, 800), label = TRUE)
#abline(v = 80, col = "grey")
abline(v = 110, col = "grey")
abline(v = 370, col = "grey")
abline(v = 430, col = "grey")
dev.off()

# Figure 5-8


pdf("xi_max.pdf", width = 5, height = 5)
mar.default <- c(5, 4, 4, 2) + 0.1
xiplot(data = data_max$Injury_Length_max, models = 30., start = 15., end = 1500., reverse = TRUE,
       ci = 0.95)
abline(h = 0, col = "grey")
abline(v = -731, col = "grey")
abline(v = -450, col = "grey")
dev.off()

# Table 5-4

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

# Table 5-5


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

# Table 5-6

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

# Figure 5-9

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

pdf("QQ_All.pdf", width = 5, height = 5)
mar.default <- c(5, 4, 4, 2) + 0.1
qqplot(x = qexp(p = ppoints(100),
                rate = 1 / beta_reg), y =  data_u$Injury_Length,
       main = "Exponential Q-Q Plot",
       xlab = "Theoretical Quantiles",
       ylab = "Sample Quantiles")
qqline(data)
dev.off()

# Figure 5-10

# This is data after subtracting the threshold 100

pdf("QQ_Max.pdf", width = 5, height = 5)
mar.default <- c(5, 4, 4, 2) + 0.1
qqplot(x = qexp(p = ppoints(100), rate = 1 / beta_reg), y = data_max_u$Injury_Length, main = "Exponential Q-Q Plot",
       xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(data)
dev.off()

