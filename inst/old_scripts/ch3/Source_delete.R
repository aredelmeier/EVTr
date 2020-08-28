
# This has *not* been updated

# This source file generates data and deletes observations that are censored
# and then calculates the MLE on the full data set (not maximums)

# generate data

gen_data <- function(censor, xi, n, num_inj, beta, rate_exp){

  set.seed(1)
  # Generate data

  Healthy <- matrix(rexp(n*num_inj, rate = rate_exp), nrow = n, ncol = num_inj)
  Injured <- matrix(QRM::rGPD(n*num_inj, xi = xi), nrow = n, ncol = num_inj)

  # Combine dataframes

  A <- matrix(NA, nrow = n, ncol = 2*num_inj)
  A[,seq(from =1, to = num_inj*2, by = 2)] <- Healthy
  A[,seq(from =2, to = num_inj*2, by = 2)] <- Injured

  # Sum injury times
  cum_sum_matrix <- t(apply(A, 1, cumsum))

  # min(cum_sum_matrix[, num_inj*2])

  B <- cum_sum_matrix - censor

  A_Ind <- ifelse(B<0,1,NA)
  B_Ind <- ifelse((B-A)*B > 0 , 0, 1)

  even <- seq(from = 2, to = 2*num_inj, 2)

  B_even <- B_Ind[,even]
  percent <- sum(B_even)/n # This gives you the percent of injuries that are censored

  temp <- A_Ind*A
  temp_injury <- temp[,even]

  # From wide to long

  Cens_long <- as.data.frame(matrix(t(temp_injury), nrow = n*num_inj, ncol = 1))
  Cens_index_long <- matrix(t(B_even), nrow = n*num_inj, ncol = 1)

  # Finalize dataframe

  Cens_long$id <- rep(1:n, each = num_inj)
  Cens_long$cens <- Cens_index_long
  Cens_long$perc <- percent
  colnames(Cens_long) <- c("Injury_Length", "ID", "Censored", "Perc_Cens")
  full <- na.omit(Cens_long)

  #    Injury_Length ID Censored Perc_Cens
  # 1      0.2211169  1        0       0.6
  # 2      1.1256608  1        0       0.6
  # 3      0.8824792  1        0       0.6
  # 4      0.4047432  1        0       0.6
  # 11     2.0664441  2        0       0.6
  # 12     0.1204254  2        0       0.6
  #
  return(full)
}


mlegpd2 <- function(data, information = c("observed", "expected"), ...) {

  injury <- data$Injury_Length
  injury <- as.numeric(injury)
  Nu <- length(injury)
  xbar <- mean(injury)
  s2 <- var(injury)

  xi0 <- -0.5 * (((xbar * xbar)/s2) - 1)
  beta0 <- 0.5 * xbar * (((xbar * xbar)/s2) + 1)
  theta <- c(xi0, beta0)

  negloglik <- function(theta, tmp) {
    xi <- theta[1]
    beta <- theta[2]
    cond1 <- beta <= 0
    cond2 <- (xi <= 0) && (max(tmp) >= (-beta/xi))
    if (cond1 || cond2)
      f <- 1e+06
    else{
      if (xi != 0){
        y <- logb(1 + (xi * tmp)/beta)
        y <- y/xi
        f <- length(tmp) * logb(beta) + (1 + xi) * sum(y)
      }
      else{
        f <- sum(logb(beta) - tmp[,1]/beta)
      }
    }
    f
  }
  fit <- optim(theta, negloglik, hessian = TRUE, ..., tmp = injury)
  return(fit$par[1])
}

mlecensgpd2 <- function (data, information = c("observed", "expected"),  ...) {

  injury <- data$Injury_Length
  cens <- data$Censored
  injury <- as.numeric(injury)

  numCens <- sum(cens)
  Nu <- length(injury)
  xbar <- mean(injury)
  s2 <- var(injury)

  xi0 <- -0.5 * (((xbar * xbar)/s2) - 1)
  beta0 <- 0.5 * xbar * (((xbar * xbar)/s2) + 1)
  theta <- c(xi0, beta0)

  data <- cbind(injury, cens)

  negloglik <- function(theta, tmp) {
    xi <- theta[1]
    beta <- theta[2]
    cond1 <- beta <= 0
    cond2 <- (xi <= 0) && (max(tmp[,1]) >= (-beta/xi))
    if (cond1 || cond2)
      f <- 1e+06
    else {
      if (xi!=0){
        y <- logb(1 + (xi * tmp[,1])/beta) # a vector
        y <- y/xi # a vector
        f <- sum((1-tmp[,2])*(logb(beta) + (1 + xi) * y) + (tmp[,2])*y)
      }
      else {
        f <- sum((1-tmp[,2])*(logb(beta) - tmp[,1]/beta) + (tmp[,2])*tmp[,1]/beta)
      }
    }
  }
  fit <- optim(theta, negloglik, ...,tmp = data)
  return(fit$par[1])
}



mle <- function(data,  method = c("MLE", "CensMLE"), information = c("observed", "expected"),  ...){
  method <- match.arg(method)
  if (method == "MLE"){
    mlegpd2(data, information = c("observed", "expected"), ...)
  }
  else{
    mlecensgpd2(data, information = c("observed", "expected"), ...)
  }
}

