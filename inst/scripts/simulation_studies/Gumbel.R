library(EVTr)

library(QRM) # use this for fit.GEV
library(evd) # use this for gevFit
library(extRemes) # use this for fevd

library(data.table)
library(ggplot2)
library(xtable)

###########################################################
############# Generate Frechet, Fit GPD ###################
###########################################################

# Truth
xi <- 0
mu <- 0
sigma <- 1
k <- 100 # 1/0.1 = 10
xk <- 5

#P(X > xk)
Prob_frechet_true <- 1 - exp(- exp(- xk))

# k-level return level
zk_frechet_true <- mu - sigma * log(- log(1 - 1 / k))

set.seed(1)
frechet <- QRM::rGEV(n = 1000, xi = xi, mu = mu, sigma = sigma)

# Fit GEV
gpd_extRemes <- extRemes::fevd(x = frechet,  type = "GEV")
gpd_extRemes_exp <- extRemes::fevd(x = frechet,  type = "Gumbel")

extRemes::lr.test(gpd_extRemes_exp, gpd_extRemes)

# Fit Gumbel!
xi_frechet_hat <- NA
xi_frechet_se <- NA
mu_frechet_hat <- as.numeric(extRemes::fevd(x = frechet,  type = "Gumbel")$results$par[1])
mu_frechet_se <- 0.03302494
sigma_frechet_hat <- as.numeric(extRemes::fevd(x = frechet,  type = "Gumbel")$results$par[2])
sigma_frechet_se <- 0.02470942

#P(X > xm)
Prob_frechet_hat <- 1 - exp(- exp(- (xk - mu_frechet_hat) / sigma_frechet_hat))

# the k-level return level is
zk_frechet_hat <- mu_frechet_hat - sigma_frechet_hat * log(- log(1 - 1 / k))

# Fit Pareto
gpd_extRemes <- extRemes::fevd(x = frechet, threshold = 0, type = "GP")
gpd_extRemes_exp <- extRemes::fevd(x = frechet, threshold = 0, type = "Exponential")

extRemes::lr.test(gpd_extRemes_exp, gpd_extRemes)

# Fit exponential!
xi_pareto_hat <- NA
xi_pareto_se <- NA
beta_pareto_hat <- as.numeric(extRemes::fevd(x = frechet, threshold = 0, type = "Exponential")$results$par)
beta_pareto_se <- 0.05016097

# P(X > xk)
Prob_pareto_hat <- exp(- xk / beta_pareto_hat)

# k-level return level
zk_pareto_hat <- beta_pareto_hat * log(k)

frechet_table <- matrix(NA, nrow = 6, ncol = 5)
rownames(frechet_table) <- c("xi", "mu", "sigma", "beta", "P(X > x_k)", "x_k")
colnames(frechet_table) <- c("Truth", "Fit_Frechet", "se", "Fit_GPD", "se")

frechet_table[, 1] <- c(xi, mu, sigma, NA, 1 / Prob_frechet_true, zk_frechet_true)
frechet_table[, 2] <- c(xi_frechet_hat, mu_frechet_hat, sigma_frechet_hat, NA, 1 / Prob_frechet_hat, zk_frechet_hat)
frechet_table[, 3] <- c(xi_frechet_se, mu_frechet_se, sigma_frechet_se, NA, NA, NA)
frechet_table[, 4] <- c(xi_pareto_hat, NA, NA, beta_pareto_hat,  1 / Prob_pareto_hat, zk_pareto_hat)
frechet_table[, 5] <- c(xi_pareto_se, NA, NA, beta_pareto_se, NA, NA)

round(frechet_table, 5)

xtable(frechet_table, digits = c(3, 3, 3, 3, 3, 3))


###########################################################
############# Generate GPD, Fit Frechet ###################
###########################################################

xi <- 0
beta <- 1
k <- 100
xk <- 5

#P(X > xk)
Prob_pareto_true <- exp(-xk / beta)

# k-level return level
zk_pareto_true <- beta * log(k)

pareto_solve <- function(beta) {

  pareto_xk <- function(xk) {
    x <- exp(-xk / beta) - 1 / k
    x
  }
  root <- uniroot(pareto_xk, lower = 0.00000001, upper = 100000)$root
  return(root)
}

set.seed(1)
pareto <- rexp(n = 1000, rate = beta)
plot(pareto)

# GEV
gpd_extRemes <- extRemes::fevd(x = pareto, type = "GEV")
gpd_extRemes_exp <- extRemes::fevd(x = pareto, type = "Gumbel")

extRemes::lr.test(gpd_extRemes_exp, gpd_extRemes)

# Fit Frechet

xi_frechet_hat2 <- as.numeric(QRM::fit.GEV(pareto)$par.ests[1])
xi_frechet_se2 <- as.numeric(QRM::fit.GEV(pareto)$par.ses[1])
mu_frechet_hat2 <- as.numeric(QRM::fit.GEV(pareto)$par.ests[2])
mu_frechet_se2 <- as.numeric(QRM::fit.GEV(pareto)$par.ses[2])
sigma_frechet_hat2 <- as.numeric(QRM::fit.GEV(pareto)$par.ests[3])
sigma_frechet_se2 <- as.numeric(QRM::fit.GEV(pareto)$par.ses[3])

# P(X > xm)
Prob_frechet_hat2 <- 1 - exp(- (1 + xi_frechet_hat2 * ((xk - mu_frechet_hat2) /
                                                         sigma_frechet_hat2)) ^ (- 1 / xi_frechet_hat2))

# the k-level return level is
zk_frechet_hat2 <- mu_frechet_hat2 - (sigma_frechet_hat2 /
                                        xi_frechet_hat2) * (1 - (-log(1 - 1 / k)) ^ (- xi_frechet_hat2))

# Fit Pareto
gpd_extRemes <- extRemes::fevd(x = pareto, threshold = 0, type = "GP")
gpd_extRemes_exp <- extRemes::fevd(x = pareto, threshold = 0, type = "Exponential")

extRemes::lr.test(gpd_extRemes_exp, gpd_extRemes)

# Fit Exponential
xi_pareto_hat2 <- NA
xi_pareto_se2 <- NA
beta_pareto_hat2 <- as.numeric(extRemes::fevd(x = pareto, threshold = 0, type = "Exponential")$results$par)
beta_pareto_se2 <- 0.03261256

# P(X > xk)
Prob_pareto_hat2 <- exp(- xk / beta_pareto_hat2)

# k-level return level
zk_pareto_hat2 <- beta_pareto_hat2 * log(k)


pareto_table <- matrix(NA, nrow = 6, ncol = 5)
rownames(pareto_table) <- c("xi", "beta", "mu", "sigma", "P(X > x_k)", "x_k")
colnames(pareto_table) <- c("Truth", "Fit_Frechet", "se", "Fit_GPD", "se")

pareto_table[, 1] <- c(xi, beta , NA, NA, 1 / Prob_pareto_true, zk_pareto_true)
pareto_table[, 2] <- c(xi_frechet_hat2, NA, mu_frechet_hat2, sigma_frechet_hat2, 1 / Prob_frechet_hat2, zk_frechet_hat2)
pareto_table[, 3] <- c(xi_frechet_se2, NA, mu_frechet_se2, sigma_frechet_se2, NA, NA)
pareto_table[, 4] <- c(xi_pareto_hat2, beta_pareto_hat2, NA, NA, 1 / Prob_pareto_hat2, zk_pareto_hat2)
pareto_table[, 5] <- c(xi_pareto_se2, beta_pareto_se2, NA, NA, NA, NA)

round(pareto_table, 3)
xtable(pareto_table, digits = c(3, 3, 3, 3, 3, 3))
