# Example

library(EVTr)
library(QRM)
library(data.table)

# parameters
n.sim <- 100
n <- 100
num_inj <- 30
xi <- 0.5
censor <- 10
rate_exp <- 1
method <- c("MLE", "CensMLE")

# 1. generate data
data <- gen_data(censor = censor,
                 xi = xi,
                 n = n,
                 num_inj = num_inj,
                 rate_exp = rate_exp,
                 specific = "keep_censored_obs",
                 seed = 10)


# regular MLE estimate
mle(data = data, method = "MLE")
# 0.1440268

# censored MLE estimate
mle(data = data, method = "CensMLE")
# 0.3094609



# 2. simulate data and fit using one function
# regular MLE
simulation(censor = censor,
           xi = xi,
           n = n,
           num_inj = num_inj,
           rate_exp = rate_exp,
           method = "MLE",
           specific = "keep_censored_obs",
           seed = 10)



# censored MLE estimate
simulation(censor = censor,
           xi = xi,
           n = n,
           num_inj = num_inj,
           rate_exp = rate_exp,
           method = "CensMLE",
           specific = "keep_censored_obs",
           seed = 10)
