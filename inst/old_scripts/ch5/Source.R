
gen_data <- function(censor, xi, n, num_inj, rate_exp){
  
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
  
  # From wide to long 
  
  Cens_long <- as.data.frame(matrix(t(temp_injury), nrow = n*num_inj, ncol = 1))
  Cens_index_long <- as.vector(t(Censored_Matrix_inj))
  
  # Finalize dataframe 
  
  Cens_long$id <- rep(1:n, each = num_inj)  
  Cens_long$cens <- Cens_index_long
  colnames(Cens_long) <- c("Injury_Length", "ID", "Censored")
  full <- na.omit(Cens_long)
  
  return(full)
}

maximum <- function(full){
  temp <- full %>% dplyr::group_by(ID) %>% dplyr::summarise(max = max(Injury_Length))
  data <- left_join(temp, full, by = c("max"="Injury_Length", "ID"="ID")) %>% dplyr::mutate(Injury_Length = max) %>% dplyr::select(-max)
  
  #df <- as.data.frame(data)
  #perc <- sum(df$Censored)/1000 # if you divide by n you get an error!
  
  #data <- data %>% dplyr::mutate(perc = perc)
  return(data)
}

excess <- function(data, ne){
  thresh <- findthreshold(data = data$Injury_Length, ne = ne)
  pot_data <- data %>% dplyr::filter(Injury_Length >=thresh) %>% dplyr::mutate(excess = Injury_Length - thresh)
  pot_data <- pot_data %>% dplyr::mutate(Injury_Length = excess, Sum = sum(Censored), perc = Sum/nrow(pot_data)) %>% dplyr::select(-Sum)
  
  return(pot_data)
}

#data <- CdS
#threshold <- 100
#information <- "observed"

#########################################################
################# MLE/CensMLE FIXED #####################
#########################################################

mlegpd <- function (data, threshold, information = c("observed", "expected"), verbose = TRUE, ...) {
  
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
  pvec <- ((1:Nu) + delta)/(Nu + delta)
  a1 <- mean(sort(excess) * (1 - pvec))
  xi <- 2 - a0/(a0 - 2 * a1)
  beta <- (2 * a0 * a1)/(a0 - 2 * a1)
  par.ests <- c(xi, beta)
  
    negloglik <- function(theta, ydata) {
      -sum(dGPD(ydata, theta[1], abs(theta[2]), log = TRUE))
    }
    
    deriv <- function(theta, ydata) {
      xi <- theta[1]
      beta <- theta[2]
      term1 <- sum(ydata/(beta + xi * ydata))
      term2 <- sum(log(1 + xi * ydata/beta))
      d1 <- -term2 * xi^(-2) + (1 + 1/xi) * term1
      d2 <- (length(ydata) - (xi + 1) * term1)/beta
      c(d1, d2)
    }
    
    fit <- optim(par.ests, fn = negloglik, gr = deriv, ydata = excess)
    
    par.ests <- fit$par
    par.ests[2] <- abs(par.ests[2])
    
    
    if (information == "observed") {
      fisher <- hessian(negloglik, fit$par, ydata = excess)
      varcov <- solve(fisher)
    }
    if (information == "expected") {
      one <- (1 + par.ests[1])^2/Nu
      two <- (2 * (1 + par.ests[1]) * par.ests[2]^2)/Nu
      cov <- -((1 + par.ests[1]) * par.ests[2])/Nu
      varcov <- matrix(c(one, cov, cov, two), 2)
    }
  
 
  
    par.ses <- sqrt(diag(varcov))
    p.less.thresh <- 1 - Nu/n
    out <- list(n = length(data), 
              p.less.thresh = p.less.thresh, n.exceed = Nu, 
              par.ests = par.ests, par.ses = par.ses, varcov = varcov 
              )
    
    names(out$par.ests) <- c("xi", "beta")
    names(out$par.ses) <- c("xi", "beta")
    out
}

mlecensgpd <- function (data, threshold, information = c("observed", "expected"), verbose = TRUE, ...) {
  
  
  information <- match.arg(information)
  
  n <- nrow(data)
  
  data <- data %>% dplyr::filter(Injury_Length > threshold)
  cens <- data$Censored
  exceedances <- data$Injury_Length
  excess <- exceedances - threshold
  Nu <- length(excess)
  xbar <- mean(excess)
  a0 <- xbar
  gamma <- -0.35
  delta <- 0
  pvec <- ((1:Nu) + delta)/(Nu + delta)
  a1 <- mean(sort(excess) * (1 - pvec))
  xi <- 2 - a0/(a0 - 2 * a1)
  beta <- (2 * a0 * a1)/(a0 - 2 * a1)
  par.ests <- c(xi, beta)

  negloglik <- function(theta, ydata) { 
    xi <- theta[1]
    beta <- theta[2]
    cond1 <- beta <= 0
    cond2 <- (xi <= 0) && (max(ydata[,1]) >= (-beta/xi))
    if (cond1 || cond2) 
      f <- 1e+06
    else { 
      y <- logb(1 + (xi * ydata[,1])/beta)
      y <- y/xi 
      f <- sum((1-ydata[,2])*(logb(beta) + (1 + xi) * y) + (ydata[,2])*y)
      
    }
  }
  
  fit <- optim(par.ests, fn = negloglik, ydata = cbind(excess,cens))
  
  par.ests <- fit$par
  par.ests[2] <- abs(par.ests[2])
  
  
  if (information == "observed") {
    fisher <- hessian(negloglik, fit$par, ydata = cbind(excess, cens))
    varcov <- solve(fisher)
  }
  if (information == "expected") {
    one <- (1 + par.ests[1])^2/Nu
    two <- (2 * (1 + par.ests[1]) * par.ests[2]^2)/Nu
    cov <- -((1 + par.ests[1]) * par.ests[2])/Nu
    varcov <- matrix(c(one, cov, cov, two), 2)
  }
  
  par.ses <- sqrt(diag(varcov))
  p.less.thresh <- 1 - Nu/n
  out <- list(n = n, 
              p.less.thresh = p.less.thresh, n.exceed = Nu, 
              par.ests = par.ests, par.ses = par.ses, varcov = varcov 
  )
  
  names(out$par.ests) <- c("xi", "beta")
  names(out$par.ses) <- c("xi", "beta")
  out
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
  fit <- optim(theta, fn = negloglik, ..., negloglik, tmp = injury)
  return(fit$par[1])
}

mlecensgpd2 <- function (data, information = c("observed", "expected"), ...) {
  
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
  fit <- optim(theta, negloglik, ..., negloglik, tmp = data)
  return(fit$par[1])
}

mle <- function(data, threshold,  method = c("MLE", "CensMLE"), information = c("observed", "expected")){
  method <- match.arg(method)
  if(method == "MLE"){
    mlegpd(data, threshold, information = c("observed", "expected"))
  }
  else {
    mlecensgpd(data, threshold, information = c("observed", "expected"))
  }

}
