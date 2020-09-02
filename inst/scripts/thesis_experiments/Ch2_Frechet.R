library(EVTr)
library(QRM) # use this for fit.GEV
library(evd) # use this for gevFit

library(data.table)
library(xtable)

# Frechet distribution = GEV distribution when xi > 0

# Gumbel distribution = GEV distribution when xi = 0

###########################################################
############# Generate Frechet, Fit GPD ###################
###########################################################

gen_Frechet_fit_GPD <- function(xi, mu, sigma, k, x_k) {

  #P(X > xk)
  Prob_frechet_true <- 1 - exp(- (1 + xi * ((xk - mu) / sigma)) ^ (- 1 / xi))

  # k-level return level
  zk_frechet_true <- mu - (sigma / xi) * (1 - (- log(1 - 1 / k)) ^ (- xi))

  set.seed(1)
  frechet <- QRM::rGEV(n = 1000, xi = xi, mu = mu, sigma = sigma)
  # plot(frechet)

  # fit with frechet
  xi_frechet_hat <- as.numeric(QRM::fit.GEV(frechet)$par.ests[1])
  xi_frechet_se <- as.numeric(QRM::fit.GEV(frechet)$par.ses[1])
  mu_frechet_hat <- as.numeric(QRM::fit.GEV(frechet)$par.ests[2])
  mu_frechet_se <- as.numeric(QRM::fit.GEV(frechet)$par.ses[2])
  sigma_frechet_hat <- as.numeric(QRM::fit.GEV(frechet)$par.ests[3])
  sigma_frechet_se <- as.numeric(QRM::fit.GEV(frechet)$par.ses[3])

  #P(X > xm)
  Prob_frechet_hat <- 1 - exp(- (1 + xi_frechet_hat * ((xk - mu_frechet_hat) /
                                                         sigma_frechet_hat)) ^ (- 1 / xi_frechet_hat))

  # the k-level return level is
  zk_frechet_hat <- mu_frechet_hat - (sigma_frechet_hat / xi_frechet_hat) * (1 - (-log(1 - 1 / k)) ^ (- xi_frechet_hat))

  # fit with Pareto
  xi_pareto_hat <- as.numeric(QRM::fit.GPD(frechet, threshold = 0)$par.ests[1])
  xi_pareto_se <- as.numeric(QRM::fit.GPD(frechet, threshold = 0)$par.ses[1])
  beta_pareto_hat <- as.numeric(QRM::fit.GPD(frechet, threshold = 0)$par.ests[2])
  beta_pareto_se <- as.numeric(QRM::fit.GPD(frechet, threshold = 0)$par.ses[2])

  #P(X > xk)
  Prob_pareto_hat <- (1 + xi_pareto_hat * xk / beta_pareto_hat) ^ (- 1 / xi_pareto_hat)

  # the k-level return level is
  zk_pareto_hat <- (beta_pareto_hat / xi_pareto_hat) * (k ^ xi_pareto_hat - 1)

  frechet_table <- matrix(NA, nrow = 6, ncol = 5)
  rownames(frechet_table) <- c("xi", "mu", "sigma", "beta", "P(X > x_k)", "x_k")
  colnames(frechet_table) <- c("Truth", "Fit_Frechet", "se", "Fit_GPD", "se")

  frechet_table[, 1] <- c(xi, mu, sigma, NA, 1 / Prob_frechet_true, zk_frechet_true)
  frechet_table[, 2] <- c(xi_frechet_hat, mu_frechet_hat, sigma_frechet_hat, NA, 1 / Prob_frechet_hat, zk_frechet_hat)
  frechet_table[, 3] <- c(xi_frechet_se, mu_frechet_se, sigma_frechet_se, NA, NA, NA)
  frechet_table[, 4] <- c(xi_pareto_hat, NA, NA, beta_pareto_hat,  1 / Prob_pareto_hat, zk_pareto_hat)
  frechet_table[, 5] <- c(xi_pareto_se, NA, NA, beta_pareto_se, NA, NA)

  # round(frechet_table, 5)
  return(frechet_table)

}

# Truth 1 - Not in Table 2-1, see Gumbel.R instead!
xi <- 0
mu <- 0
sigma <- 1
k <- 100
xk <- 5

frechet_table <- gen_Frechet_fit_GPD(xi, mu, sigma, k, x_k)
round(frechet_table, 5)
xtable::xtable(frechet_table, digits = c(3, 3, 3, 3, 3, 3))


# Truth 2 - Table 2-1, GEV RVs Generated 2nd part
xi <- 0.2
mu <- 0
sigma <- 1
k <- 100
xk <- 5

frechet_table <- gen_Frechet_fit_GPD(xi, mu, sigma, k, x_k)
round(frechet_table, 5)
xtable::xtable(frechet_table, digits = c(3, 3, 3, 3, 3, 3))



# Truth 3  - Table 2-1, GEV RVs Generated 3rd part
xi <- 0.5
mu <- 0
sigma <- 1
k <- 100
xk <- 5

###########################################################
############# Generate GPD, Fit Frechet ###################
###########################################################

gen_GPD_fit_Frechet <- function(xi, beta, k, x_k) {

  # P(X > xk)
  Prob_pareto_true <- (1 + xi * xk / beta) ^ (- 1 / xi)

  # k-level return level
  zk_pareto_true <- (beta / xi) * (k ^ xi - 1)


  set.seed(1)
  pareto <- QRM::rGPD(n = 1000, xi = xi, beta = beta)
  plot(pareto)

  # fit with frechet
  xi_frechet_hat2 <- as.numeric(QRM::fit.GEV(pareto)$par.ests[1])
  xi_frechet_se2 <- as.numeric(QRM::fit.GEV(pareto)$par.ses[1])
  mu_frechet_hat2 <- as.numeric(QRM::fit.GEV(pareto)$par.ests[2])
  mu_frechet_se2 <- as.numeric(QRM::fit.GEV(pareto)$par.ses[2])
  sigma_frechet_hat2 <- as.numeric(QRM::fit.GEV(pareto)$par.ests[3])
  sigma_frechet_se2 <- as.numeric(QRM::fit.GEV(pareto)$par.ses[3])

  #P(X > xm)
  Prob_frechet_hat2 <- 1 - exp(- (1 + xi_frechet_hat2 * ((xk - mu_frechet_hat2) /
                                                           sigma_frechet_hat2)) ^ (- 1 / xi_frechet_hat2))

  # the k-level return level is
  zk_frechet_hat2 <- mu_frechet_hat2 - (sigma_frechet_hat2 /
                                          xi_frechet_hat2) * (1 - (-log(1 - 1 / k)) ^ (- xi_frechet_hat2))

  # fit with Pareto
  xi_pareto_hat2 <- as.numeric(QRM::fit.GPD(pareto, threshold = 0)$par.ests[1])
  xi_pareto_se2 <- as.numeric(QRM::fit.GPD(pareto, threshold = 0)$par.ses[1])
  beta_pareto_hat2 <- as.numeric(QRM::fit.GPD(pareto, threshold = 0)$par.ests[2])
  beta_pareto_se2 <- as.numeric(QRM::fit.GPD(pareto, threshold = 0)$par.ses[2])

  #P(X > xk)
  Prob_pareto_hat2 <- (1 + xi_pareto_hat2 * xk / beta_pareto_hat2) ^ (- 1 / xi_pareto_hat2)

  # the k-level return level is
  zk_pareto_hat2 <- (beta_pareto_hat2 / xi_pareto_hat2) * (k ^ xi_pareto_hat2 - 1)

  pareto_table <- matrix(NA, nrow = 6, ncol = 5)
  rownames(pareto_table) <- c("xi", "beta", "mu", "sigma", "P(X > x_k)", "x_k")
  colnames(pareto_table) <- c("Truth", "Fit_Frechet", "se", "Fit_GPD", "se")

  pareto_table[, 1] <- c(xi, beta, NA, NA, 1 / Prob_pareto_true, zk_pareto_true)
  pareto_table[, 2] <- c(xi_frechet_hat2, NA, mu_frechet_hat2, sigma_frechet_hat2,
                         1 / Prob_frechet_hat2, zk_frechet_hat2)
  pareto_table[, 3] <- c(xi_frechet_se2, NA, mu_frechet_se2, sigma_frechet_se2, NA, NA)
  pareto_table[, 4] <- c(xi_pareto_hat2, beta_pareto_hat2, NA, NA, 1 / Prob_pareto_hat2,
                         zk_pareto_hat2)
  pareto_table[, 5] <- c(xi_pareto_se2, beta_pareto_se2, NA, NA, NA, NA)

  # round(pareto_table, 3)
  return(pareto_table)
}

# Truth 1 - Not in Table 2-1, see Gumbel.R instead!
xi <- 0
beta <- 1
k <- 100
xk <- 5

pareto_table <- gen_GPD_fit_frechet_and_GP(xi, beta, k, xk)

round(pareto_table, 3)

xtable(pareto_table, digits = c(3, 3, 3, 3, 3, 3))

# Truth 2 - Table 2-1, GP RVs Generated 2nd part
xi <- 0.2
beta <- 1
k <- 100
xk <- 5

pareto_table <- gen_GPD_fit_Frechet(xi, beta, k, xk)

round(pareto_table, 3)

xtable(pareto_table, digits = c(3, 3, 3, 3, 3, 3))


# Truth 2 - Table 2-1, GP RVs Generated 3rd part
xi <- 0.5
beta <- 1
k <- 100
xk <- 5

pareto_table <- gen_GPD_fit_Frechet(xi, beta, k, xk)

round(pareto_table, 3)

xtable(pareto_table, digits = c(3, 3, 3, 3, 3, 3))
