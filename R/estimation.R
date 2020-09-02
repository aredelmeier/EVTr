#' A function to estimate the maximum likelihood estimator for either the regular generalized Pareto
#' distribution or the censored generalized Pareto distribution.
#'
#' @param data Dataframe. One column is called "Injury_Length" and corresponds to the observations.
#' If method is either "CensMLE" or "CensMLE_full", the second column must be a binary indicator of
#' whether the injuries (the first column) are censored or not.
#'
#' @param method String. One of "MLE", "CensMLE", "MLE_full" or "CensMLE_full".
#'
#' @param threshold Numeric. (Optional). Used only when method = "MLE_full" or "CensMLE_full".
#' This determines the exceedances to use when estimating the parameters of the generalzied Pareto
#' distribution.
#'
#' @param information String. (Optional). Used only when method = "MLE_full" or "CensMLE_full".
#' This determines the information criterion to use when optimizing.
#'
#' @return Numeric. Corresponding to the maximum likelihood estimator based on one of the generalized
#' Pareto distribution.
#'
#' @export
#'
#' @author Annabelle Redelmeier
#'
#' @examples
#'
#' data <- data.frame(Injury_Length = rexp(10), Censored = rbinom(10, 1, 0.5))
#'
#' mle(data = data, method = "CensMLE")
#'
#' mle(data = data, method = "CensMLE", threshold = 0, information = "observed")
mle <- function(data,
                method = c("MLE",
                           "CensMLE",
                           "MLE_full",
                           "CensMLE_full"),
                threshold = NULL,
                information = c("observed", "expected")) {

  method <- match.arg(method)
  if (method == "MLE") {
    mle_gpd(data)

  } else if (method == "CensMLE") {
    mle_cens_gpd(data)

  } else if (method == "MLE_full") {

    information <- match.arg(information)
    mle_gpd_full(data = data, threshold = threshold, information = information)

  } else if (method == "CensMLE_full") {

    information <- match.arg(information)
    mle_cens_gpd_full(data = data, threshold = threshold, information = information)
  }


}

#' A function that calculates the maximum likelihood estimator based on the generalized Pareto distribution.
#'
#' @param data Dataframe. One column is called "Injury_Length" and corresponds to the observations.

#' @return Numeric. Corresponding to the maximum likelihood estimator based on the generalized Pareto
#' distribution.
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

  xi0 <- -0.5 * (((xbar * xbar) / s2) - 1)
  beta0 <- 0.5 * xbar * (((xbar * xbar) / s2) + 1)
  theta <- c(xi0, beta0)


  fit <- optim(theta, negloglik_gpd, method = "BFGS", hessian = FALSE, tmp = injury)
  return(fit$par[1])
}

#' The negative log-likelihood of the generalized Pareto distribution.
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
  cond2 <- (xi <= 0) && (max(tmp) >= (- beta / xi))
  if (cond1 || cond2)
    f <- 1e+06
  else{
    if (xi != 0) {
      y <- logb(1 + (xi * tmp) / beta)
      y <- y / xi
      f <- length(tmp) * logb(beta) + (1 + xi) * sum(y)
    }
    else{
      f <- sum(logb(beta) - tmp[, 1] / beta)
    }
  }
  f
}

#' A function that calculates the maximum likelihood estimator based on the censored generalized Pareto
#' distribution.
#'
#' @param data Dataframe. One column is called "Injury_Length" and corresponds to the observations. One
#' column called "Censored" of 0s and 1s. 0 indicates that the observation is not censored. 1 indicates the
#' observation is censored.
#'
#' @return Numeric. Corresponding to the maximum likelihood estimator based on the censored generalized
#' Pareto distribution.
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

  xi0 <- -0.5 * (((xbar * xbar) / s2) - 1)
  beta0 <- 0.5 * xbar * (((xbar * xbar) / s2) + 1)
  theta <- c(xi0, beta0)

  data <- cbind(injury, cens)

  fit <- optim(theta, negloglik_cens_gpd, method = "BFGS", hessian = FALSE, tmp = data)
  return(fit$par[1])
}

#' The negative censored log-likelihood of the generalized Pareto distribution.
#'
#' @param theta Vector. Length 2. The first argument corresponds to the xi parameter in the GPD distribution and
#' the second argument corresponds to the beta parameter.
#'
#' @param tmp Dataframe. Two columns. First column are the observations. The second column is binary indicating
#' which of the observations are censored.
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
  cond2 <- (xi <= 0) && (max(tmp[, 1]) >= (- beta / xi))
  if (cond1 || cond2)
    f <- 1e+06
  else {
    if (xi != 0) {
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


#' This function actually comes from the R package QRM and is only used to test to see if the other mle
#' functions work properly.
#'
#' @param data Dataframe. One column is called "Injury_Length" and corresponds to the observations.
#'
#' @param threshold Numeric. Threshold to use when using the peaks-over-threshold method.
#'
#' @param information String. One of "observed" or "expected" used in optimization.
#'
#' @return List. Contains  "n" which is the number of observations in \code{data},
#' "p.less.thresh" which is the proportion below the threshold,
#' "n.exceed" is the number of observations above the threshold, "par.ests" is a vector
#' containing the parameter estimates of xi and beta, "par.ses" contains the standard errors
#' of the estimated parameters, and "varcov" contains the covariance matrix of the estimated
#' parameters.
#'
#' @export
#'
#' @import numDeriv
#'
#' @examples
#'
#' data <- data.frame(Injury_Length = rexp(100))
#'
#' mle_gpd_full(data = data, threshold = 1, information = "expected")
mle_gpd_full <- function(data, threshold, information = c("observed", "expected")) {

  data <- data$Injury_Length

  information <- match.arg(information)

  n <- length(data)

  exceedances <- data[data > threshold]
  excess <- exceedances - threshold
  Nu <- length(excess)
  xbar <- mean(excess)
  a0 <- xbar
  gamma <- -0.35
  delta <- 0
  pvec <- ((1:Nu) + delta) / (Nu + delta)
  a1 <- mean(sort(excess) * (1 - pvec))
  xi <- 2 - a0 / (a0 - 2 * a1)
  beta <- (2 * a0 * a1) / (a0 - 2 * a1)
  par.ests <- c(xi, beta)

  negloglik <- function(theta, ydata) {
    -sum(dGPD_QRM(ydata, theta[1], abs(theta[2]), log = TRUE))
  }

  deriv <- function(theta, ydata) {
    xi <- theta[1]
    beta <- theta[2]
    term1 <- sum(ydata / (beta + xi * ydata))
    term2 <- sum(log(1 + xi * ydata / beta))
    d1 <- - term2 * xi ^ (-2) + (1 + 1 / xi) * term1
    d2 <- (length(ydata) - (xi + 1) * term1) / beta
    return(c(d1, d2))
  }

  fit <- optim(par.ests, fn = negloglik, gr = deriv, ydata = excess)

  par.ests <- fit$par
  par.ests[2] <- abs(par.ests[2])


  if (information == "observed") {
    fisher <- hessian(negloglik, fit$par, ydata = excess)
    varcov <- solve(fisher)
  }
  if (information == "expected") {
    one <- (1 + par.ests[1]) ^ 2 / Nu
    two <- (2 * (1 + par.ests[1]) * par.ests[2] ^ 2) / Nu
    cov <- - ((1 + par.ests[1]) * par.ests[2]) / Nu
    varcov <- matrix(c(one, cov, cov, two), 2)
  }

  par.ses <- sqrt(diag(varcov))
  p.less.thresh <- 1 - Nu / n
  out <- list(n = length(data),
              p.less.thresh = p.less.thresh, n.exceed = Nu,
              par.ests = par.ests, par.ses = par.ses, varcov = varcov
  )

  names(out$par.ests) <- c("xi", "beta")
  names(out$par.ses) <- c("xi", "beta")
  out
}

#' This function is the dGPD function from the QRM package.
#'
#' @param x Vector. The data used to estimate GPD parameters from.
#'
#' @param xi Numeric. The xi parameter from the GPD distribution.
#'
#' @param beta Numeric. The beta parameter from the GPD distribution.
#'
#' @param log Logical. If TRUE, will not return the log of the estimate.
#'
#' @return Vector of length \code{length(x)}. The density of the GPD evaluated at
#' each \code{x} value.
#'
#' @export
#'
#' @examples
#'
#' x <- 1:10; xi <- 0.5;
#'
#' dGPD_QRM(x = x, xi = xi)
dGPD_QRM <- function(x, xi, beta = 1, log = FALSE) {
  xb <- x / beta

  if (xi == 0) {
    res <- log(dexp(xb)) - log(beta)
  } else {

    if (xi < 0) {
      ind <- xb > 0 & (xb < 1 / abs(xi))
    } else {
      ind <- xb > 0
    }

    r <- rep(-Inf, length(x))
    r[ind] <- (-1 / xi - 1) * log(1 + xi * xb[ind]) - log(beta)
    res <- r
  }

  if (log) {
    return(res)
  } else {
    return(exp(res))
  }
}

#' This function is a mixture of two packages (I think): one that does mle estimation and the other
#' that defines the censored negative log likelihood of the generalized Pareto distribution.
#'
#' @param data Dataframe. One column is called "Injury_Length" and corresponds to the observations. One
#' column called "Censored" of 0s and 1s. 0 indicates that the observation is not censored. 1 indicates the
#' observation is censored.
#'
#' @param threshold Numeric. Threshold to use when using the peaks-over-threshold method.
#'
#' @param information String. One of "observed" or "expected" used in optimization.
#'
#' @return List. Contains  "n" which is the number of observations in \code{data},
#' "p.less.thresh" which is the proportion below the threshold,
#' "n.exceed" is the number of observations above the threshold, "par.ests" is a vector
#' containing the parameter estimates of xi and beta, "par.ses" contains the standard errors
#' of the estimated parameters, and "varcov" contains the covariance matrix of the estimated
#' parameters.
#'
#' @export
#'
#' @author Annabelle Redelmeier
#'
#' @examples
#'
#' data <- data.frame(Injury_Length = rexp(100), Censored = rbinom(100, 1, 0.5))
#'
#' mle_cens_gpd_full(data = data, threshold = 1, information = "expected")
mle_cens_gpd_full <- function(data, threshold, information = c("observed", "expected")) {

  Injury_Length <- NULL

  data <- data.table(data)

  information <- match.arg(information)

  n <- nrow(data)

  data <- data[Injury_Length > threshold]
  cens <- data$Censored
  exceedances <- data$Injury_Length
  excess <- exceedances - threshold
  Nu <- length(excess)
  xbar <- mean(excess)
  a0 <- xbar
  gamma <- -0.35
  delta <- 0
  pvec <- ((1:Nu) + delta) / (Nu + delta)
  a1 <- mean(sort(excess) * (1 - pvec))
  xi <- 2 - a0 / (a0 - 2 * a1)
  beta <- (2 * a0 * a1) / (a0 - 2 * a1)
  par.ests <- c(xi, beta)

  negloglik <- function(theta, ydata) {
    xi <- theta[1]
    beta <- theta[2]
    cond1 <- beta <= 0
    cond2 <- (xi <= 0) && (max(ydata[, 1]) >= (- beta / xi))
    if (cond1 || cond2)
      f <- 1e+06
    else {
      y <- logb(1 + (xi * ydata[, 1]) / beta)
      y <- y / xi
      f <- sum((1 - ydata[, 2]) * (logb(beta) + (1 + xi) * y) + (ydata[, 2]) * y)

    }
  }

  fit <- optim(par.ests, fn = negloglik, ydata = cbind(excess, cens))

  par.ests <- fit$par
  par.ests[2] <- abs(par.ests[2])


  if (information == "observed") {
    fisher <- hessian(negloglik, fit$par, ydata = cbind(excess, cens))
    varcov <- solve(fisher)
  }
  if (information == "expected") {
    one <- (1 + par.ests[1]) ^ 2 / Nu
    two <- (2 * (1 + par.ests[1]) * par.ests[2] ^ 2) / Nu
    cov <- - ((1 + par.ests[1]) * par.ests[2]) / Nu
    varcov <- matrix(c(one, cov, cov, two), 2)
  }

  par.ses <- sqrt(diag(varcov))
  p.less.thresh <- 1 - Nu / n
  out <- list(n = n,
              p.less.thresh = p.less.thresh, n.exceed = Nu,
              par.ests = par.ests, par.ses = par.ses, varcov = varcov
  )

  names(out$par.ests) <- c("xi", "beta")
  names(out$par.ses) <- c("xi", "beta")
  out
}
