
#' The negative log-likelihood of the generalized Pareto distribution. Based on dGPD from the QRM R package.
#'
#' @param theta Vector. Length 2. The first argument corresponds to $$xi$$ in the GPD distribution and the second
#' argument corresponds to $$beta$$.
#'
#' @param tmp Dataframe. One column corresponding to the data/observations.
#'
#' @return Numeric. Corresponding to the maximum likelihood estimator based on the generalized Pareto distribution.
#'
#' @export
#'
#' @import stats
#'
#' @author Annabelle Redelmeier
#'
#' @examples
#'
#' data <- rexp(10)
#'
#' negloglik_gpd(theta = c(1, 1), tmp = data)
negloglik_gpd <- function(theta, tmp) {
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

#' The negative censored log-likelihood of the generalized Pareto distribution.
#'
#' @param theta Vector. Length 2. The first argument corresponds to $$xi$$ in the GPD distribution and the second
#' argument corresponds to $$beta$$.
#'
#' @param tmp Dataframe. Two columns. First column is the data/observations. The second column is binary where 1
#' means the observation is censored and 0 means the observation is not censored.
#'
#' @return Numeric. Corresponding to the maximum likelihood estimator based on the censored generalized Pareto
#' distribution.
#'
#' @export
#'
#' @author Annabelle Redelmeier
#'
#' @examples
#'
#' data <- rexp(10)
#' data <- cbind(data, censored = rbinom(10, 1, 0.5))
#'
#' negloglik_cens_gpd(theta = c(1, 1), tmp = data)
negloglik_cens_gpd <- function(theta, tmp) {
  xi <- theta[1]
  beta <- theta[2]
  cond1 <- beta <= 0
  cond2 <- (xi <= 0) && (max(tmp[, 1]) >= (-beta / xi))
  if (cond1 || cond2)
    f <- 1e+06
  else {
    if (xi!=0){
      y <- logb(1 + (xi * tmp[, 1]) / beta)
      y <- y / xi
      f <- sum((1 - tmp[, 2]) * (logb(beta) + (1 + xi) * y) + (tmp[, 2]) * y)
    }
    else {
      f <- sum((1 - tmp[, 2]) * (logb(beta) - tmp[, 1] / beta) + (tmp[, 2]) * tmp[, 1] / beta)
    }
  }
  f
}

#' The negative log-likelihood of the generalized Pareto distribution implemented in the \code{QRM} R package.
#'
#' @param theta Vector. Length 2. The first argument corresponds to $$xi$$ in the GPD distribution and the second
#' argument corresponds to $$beta$$.
#'
#' @param ydata Dataframe. One column corresponding to the data/observations.
#'
#' @return Numeric. Corresponding to the maximum likelihood estimator based on the generalized Pareto distribution.
#'
#' @export
#'
#' @author Annabelle Redelmeier
#'
#' @examples
#'
#' data <- rexp(10)
#'
#' negloglik_gpd_full(theta = c(1, 1), ydata = data)
negloglik_gpd_full <- function(theta, ydata) {
  -sum(dGPD(ydata, theta[1], abs(theta[2]), log = TRUE))
}


#' The negative censored log-likelihood of the censored generalized Pareto distribution implemented by the \code{dGPD} R package.
#'
#' @param theta Vector. Length 2. The first argument corresponds to $$xi$$ in the GPD distribution and the second
#' argument corresponds to $$beta$$.
#'
#' @param ydata Dataframe. the first column corresponding to the data/observations. The second column corresponding to the
#' binary column. 0 if the observation is not censored. 1 if the observation is censored.
#'
#' @return Numeric. Corresponding to the maximum likelihood estimator based on the generalized Pareto distribution.
#'
#' @export
#'
#' @author Annabelle Redelmeier
#'
#' @examples
#'
#' data <- data.frame(obs = rexp(10), Censored = rbinom(10, 1, 0.5))
#'
#' negloglik_cens_gpd_full(theta = c(1, 1), ydata = data)

negloglik_cens_gpd_full <- function(theta, ydata) {
  xi <- theta[1]
  beta <- theta[2]
  cond1 <- beta <= 0
  cond2 <- (xi <= 0) && (max(ydata[,1]) >= (-beta/xi))
  if (cond1 || cond2)
    f <- 1e+06
  else {
    y <- logb(1 + (xi * ydata[,1])/beta)
    y <- y/xi
    f <- sum((1 - ydata[, 2]) * (logb(beta) + (1 + xi) * y) + (ydata[, 2]) * y)

  }
  return(f)
}

#' A function that calculates the maximum likelihood estimator based on the generalized Pareto distribution.
#'
#' @param data Dataframe. One column is called "Injury_Length" and corresponds to the observations.

#' @return Numeric. Corresponding to the maximum likelihood estimator based on the generalized Pareto distribution.
#'
#' @export
#'
#' @author Annabelle Redelmeier
#'
#' @examples
#'
#' data <- data.frame(Injury_Length = rexp(10))
#'
#' mle_gpd(data = data)
mle_gpd <- function(data) {

  var <- optim <- NULL

  injury <- data$Injury_Length
  injury <- as.numeric(injury)
  Nu <- length(injury)
  xbar <- mean(injury)
  s2 <- var(injury)

  xi0 <- -0.5 * (((xbar * xbar)/s2) - 1)
  beta0 <- 0.5 * xbar * (((xbar * xbar)/s2) + 1)
  theta <- c(xi0, beta0)


  fit <- optim(theta, negloglik_gpd, hessian = TRUE, tmp = injury)
  return(fit$par[1])
}

#' A function that calculates the maximum likelihood estimator based on the censored generalized Pareto distribution.
#'
#' @param data Dataframe. One column is called "Injury_Length" and corresponds to the observations. One column called
#' "Censored" of 0s and 1s. 0 indicates that the observation is not censored. 1 indicates the observation is censored.

#' @return Numeric. Corresponding to the maximum likelihood estimator based on the censored generalized Pareto distribution.
#'
#' @export
#'
#' @author Annabelle Redelmeier
#'
#' @examples
#'
#' data <- data.frame(Injury_Length = rexp(10), Censored = rbinom(10, 1, 0.5))
#'
#' mle_cens_gpd(data = data)
mle_cens_gpd <- function(data) {

  var <- optim <- NULL

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


  fit <- optim(theta, negloglik_cens_gpd, tmp = data)
  return(fit$par[1])
}

#' A function to estimate the maximum likelihood estimator for a variety of methods.
#'
#' @param data Dataframe. One column is called "Injury_Length" and corresponds to the observations.
#'
#' @param method String. One of "MLE", "CensMLE", "MLE_full", or "CensMLE_full".
#'
#' @return Numeric. Corresponding to the maximum likelihood estimator based on one of the generalized Pareto distribution.
#'
#' @export
#'
#' @author Annabelle Redelmeier
#'
#' @examples
#'
#' data <- data.frame(Injury_Length = rexp(10), Censored = rbinom(10, 1, 0.5))
#'
#' mle(data = data, method = "MLE")

mle <- function(data,
                method = c("MLE",
                           "CensMLE")){

  method <- match.arg(method)
  if (method == "MLE") {

    mle_gpd(data)

  } else if (method == "CensMLE") {

    mle_cens_gpd(data)
  }
}
