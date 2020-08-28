
# This source file generates data and keeps censored observations (part that is not censored)
# and then calculates the MLE on the maximums

gen_data <- function(censor, xi, n, num_inj, beta, rate_exp){

  set.seed(1)

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

  # Compute actual length of censored injury
  D_All <- A*(A_Ind + B_Ind)
  D_temp <- ifelse(D_All == 0, NA, D_All)
  D_even <- D_temp[,even]
  D_long <- Cens_index_long <- as.vector(t(D_even))

  # Compute indicator matrix of censored times
  C_temp <- ifelse(C == 0, 0, NA)
  Censored_Matrix <- B_Ind + C_temp
  Censored_Matrix_inj <- Censored_Matrix[,even]

  # Need a way to put a 1 if any injury censored 0 otherwise
  E_temp <- B_Ind[,even]
  id_matrix <- matrix(1, nrow = num_inj, ncol = 1)
  Cens_Ind <- E_temp%*%id_matrix
  Cens_Ind_All <- rep(Cens_Ind, each = num_inj)

  # From wide to long
  Cens_long <- as.data.frame(matrix(t(temp_injury), nrow = n*num_inj, ncol = 1))
  Cens_index_long <- as.vector(t(Censored_Matrix_inj))

  # Finalize dataframe
  Cens_long$id <- rep(1:n, each = num_inj)
  Cens_long$cens <- Cens_index_long
  Cens_long$actual <- D_long
  Cens_long$Sum_Censored <- Cens_Ind_All

  colnames(Cens_long) <- c("Injury_Length", "ID", "Censored", "Actual", "Any_Injury_Censored")
  full <- na.omit(Cens_long)

  #    Injury_Length ID Censored    Actual Any_Injury_Censored
  # 1      0.2211169  1        0 0.2211169                   1
  # 2      1.1256608  1        0 1.1256608                   1
  # 3      0.8824792  1        0 0.8824792                   1
  # 4      0.4047432  1        0 0.4047432                   1
  # 5      0.3404012  1        1 0.8599449                   1
  # 11     2.0664441  2        0 2.0664441                   0

  return(full)
}


maximum <- function(full){
  temp <- full %>% dplyr::group_by(ID) %>% dplyr::summarise(max = max(Injury_Length))
  data <- left_join(temp, full, by = c("max"="Injury_Length", "ID"="ID")) %>% dplyr::mutate(Injury_Length = max) %>% dplyr::select(-max)

  df <- as.data.frame(data)
  perc <- sum(df$Censored)/1000 # if you divide by n you get an error!

  data <- data %>% dplyr::mutate(perc = perc, Not_Censored = 1 - Censored, delta = 1-Not_Censored*Any_Injury_Censored) %>%
    dplyr::select(ID, Injury_Length, Actual,Censored, Not_Censored,Any_Injury_Censored, delta)

  return(data)
}

excess <- function(data, ne){
  thresh <- findthreshold(data = data$Injury_Length, ne = ne)
  pot_data <- data %>% dplyr::filter(Injury_Length >=thresh) %>% dplyr::mutate(excess = Injury_Length - thresh)
  pot_data <- pot_data %>% dplyr::mutate(Injury_Length_before = Injury_Length, Injury_Length = excess) %>% dplyr::select(-excess)
  return(pot_data)
}

# N calculates the number of not censored observations
# Norm_Value calculates the mean of the not censored injuries (can't use this for the maxima since this is the excess!)
# Cens_Value calculates the mean of the censored injuries
# Norm_Value_max calculates th mean of the not censored injuries (before subtracting the threshold)


mle <- function(data,  method = c("N", "Norm_Value", "Cens_Value", "Norm_Value_max")){
  method <- match.arg(method)

  if(method == "N"){
    nc <- data %>% dplyr::filter(Censored == 0)
    nrow(nc)
  }
  else if(method == "Norm_Value"){
    nc <- data %>% dplyr::filter(Censored == 0)
    mean(nc$Injury_Length)
  }
  else if(method == "Cens_Value"){
    C <- data %>% dplyr::filter(Censored == 1)
    mean(C$Actual)
  }
  else{
    BC <- data %>% dplyr::filter(Censored ==0)
    mean(BC$Injury_Length_before)
  }
}




