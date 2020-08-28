
# This has been updated

# This source file generates data and keeps censored observations (part that is not censored)
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

  B <- cum_sum_matrix - censor
  A_Ind <- ifelse(B<0,1,0)
  B_Ind <- ifelse((B-A)*B > 0 , 0, 1)

  C_Ind <- 1 - (A_Ind + B_Ind)
  C <- ifelse(C_Ind == 0, 0, NA)

  # Need to calculate the censored length
  # cum_sum - A will give you the end of the last injury so censor - this will give you the size until end of study

  back <- censor - (cum_sum_matrix - A)

  temp <- A*A_Ind + back*B_Ind + C

  even <- seq(from = 2, to = 2*num_inj, 2)
  temp_injury <- temp[,even]

  # Compute indicator matrix of censored times
  C_temp <- ifelse(C == 0, 0, NA)
  Censored_Matrix <- B_Ind + C_temp
  Censored_Matrix_inj <- Censored_Matrix[,even]

  B_even <- B_Ind[,even]
  #percent <- sum(B_even)/n # This gives you the percent of LAST injuries that are censored i.e # cens / 1000


  # From wide to long

  Cens_long <- as.data.frame(matrix(t(temp_injury), nrow = n*num_inj, ncol = 1))
  Cens_index_long <- matrix(t(Censored_Matrix_inj), nrow = n*num_inj, ncol = 1)



  # Finalize dataframe

  Cens_long$id <- rep(1:n, each = num_inj)
  Cens_long$cens <- Cens_index_long
  Cens_long$perc <- 1
  colnames(Cens_long) <- c("Injury_Length", "ID", "Censored", "Perc_Cens")
  full <- na.omit(Cens_long)
  percent <- sum(B_even)/nrow(full) # This gives you the percent of injuries that are censored  i.e # cens injuries / total num of injuries
  full$percent <- percent

  #    Injury_Length ID Censored Perc_Cens percent
  # 1      0.2211169  1        0         1    0.15
  # 2      1.1256608  1        0         1    0.15
  # 3      0.8824792  1        0         1    0.15
  # 4      0.4047432  1        0         1    0.15
  # 5      0.3404012  1        1         1    0.15 # CORRECT
  # 11     2.0664441  2        0         1    0.15

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

  # theta is a two parameter vector
  # tmp is the data

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

  # theta is a two parameter vector
  # tmp is the data and the censoring indicator

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



