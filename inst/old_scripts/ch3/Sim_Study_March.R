setwd("/Users/Redelmeier/Dropbox/Research/Simulations/Feb 6 New Code")

library(simsalapar)
library(parallel)
library(foreach)
library(doParallel)
library(QRM)
library(dplyr)
library(xtable)


###############################################################
####### Study 1: without censoring  ###########################
###############################################################

setwd("/Users/Redelmeier/Dropbox/Research/Simulations/Feb 6 New Code")

source("/Users/Redelmeier/Dropbox/Research/Simulations/Feb 6 New Code/Source_keep_max.R")

varList <- 
  varlist(
    n.sim = list(type = "N", expr = quote(N[sim]), value = 1000),
    n = list(type = "frozen",  value = 1000),
    num_inj = list(type = "frozen", value = 30),
    xi = list(type = "grid", expr = quote(xi),   value = c(0.1, 0.5, 1.5)),
    rate_exp = list(type = "frozen",  value = 1),
    ne = list(type = "frozen",  value = 50),
    method = list(type = "inner", expr = quote(method), value = c("MLE", "CensMLE")) 
  )


doOne <- function(n, xi, rate_exp, num_inj, method, ne){
  stopifnot(require(QRM), require(dplyr), sapply(list(n, xi, rate_exp, num_inj, ne), is.numeric))
  met <- rep(NA, length(method))
  names(met) <- method
  
  full <- gen_data_nocens(xi, n, num_inj, rate_exp)
  
  m <- maximum(full=full)
  
  exc <- excess(data = m, ne = ne)
  
  for(i in 1:length(method)){
    met[i] <-  mle(exc, method = method[i])
  }
  met
}

# Parallel computing
no_cores <- detectCores()-1
cl <- makeCluster(no_cores)
clusterExport(cl,"df")
clusterEvalQ(cl,library(QRM))
clusterEvalQ(cl,source("/Users/Redelmeier/Dropbox/Research/Simulations/Feb 6 New Code/Source_keep_max.R"))

res <- doForeach(varList, cluster=cl, seed="seq", doOne=doOne,  sfile="max_no_cens.rds", monitor=TRUE)
val <- getArray(res)
df <- array2df(val)

df_MLE_1 <- df %>% dplyr::filter(method == "MLE") %>% dplyr::filter(xi == 0.1) %>% dplyr::summarise(bias = mean(value) - 0.1, sd = sd(value)) %>% 
  dplyr::mutate(MSE = bias^2 + sd^2)
df_MLE_5 <- df %>% dplyr::filter(method == "MLE") %>% dplyr::filter(xi == 0.5) %>% dplyr::summarise(bias = mean(value) - 0.5, sd = sd(value))%>% 
  dplyr::mutate(MSE = bias^2 + sd^2)
df_MLE_15 <- df %>% dplyr::filter(method == "MLE") %>% dplyr::filter(xi == 1.5) %>% dplyr::summarise(bias = mean(value) - 1.5, sd = sd(value))%>% 
  dplyr::mutate(MSE = bias^2 + sd^2)
df_Cens_1 <- df %>% dplyr::filter(method == "CensMLE") %>% dplyr::filter(xi == 0.1) %>% dplyr::summarise(bias = mean(value) - 0.1, sd = sd(value))%>% 
  dplyr::mutate(MSE = bias^2 + sd^2)
df_Cens_5 <- df %>% dplyr::filter(method == "CensMLE") %>% dplyr::filter(xi == 0.5) %>% dplyr::summarise(bias = mean(value) - 0.5, sd = sd(value))%>% 
  dplyr::mutate(MSE = bias^2 + sd^2)
df_Cens_15 <- df %>% dplyr::filter(method == "CensMLE") %>% dplyr::filter(xi == 1.5) %>% dplyr::summarise(bias = mean(value) - 1.5, sd = sd(value))%>% 
  dplyr::mutate(MSE = bias^2 + sd^2)

MSE_NoCens <- matrix(NA, nrow = 3, ncol = 4)
colnames(MSE_NoCens) <- c("xi", "Bias", "sd", "MSE")
MSE_NoCens[,1] <- c(0.1, 0.5, 1.5)
MSE_NoCens[1,2] <- df_MLE_1[[1]]
MSE_NoCens[2,2] <- df_MLE_5[[1]]
MSE_NoCens[3,2] <- df_MLE_15[[1]]
MSE_NoCens[1,3] <- df_MLE_1[[2]]
MSE_NoCens[2,3] <- df_MLE_5[[2]]
MSE_NoCens[3,3] <- df_MLE_15[[2]]
MSE_NoCens[1,4] <- df_MLE_1[[3]]
MSE_NoCens[2,4] <- df_MLE_5[[3]]
MSE_NoCens[3,4] <- df_MLE_15[[3]]


xtable(MSE_NoCens, digits = c(1, 1, 3, 3, 3))


#setwd("/Users/Redelmeier/Dropbox/Research/Presentations/To_Show_Plots")

setwd("/Users/Redelmeier/Dropbox/Research/Thesis")

pdf("Sim_March_No_Cens.pdf",width=5, height=5)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 5, 0, 0))
boxplot(I(value-as.numeric(as.character(df$xi)))~xi, data = df)
abline(h=0,col="red")
mtext(text="True xi", side=1, line=4, las=1)
mtext(text="Observed Bias", side=2, line=4, las=0)
dev.off()


###############################################################
####### Study 1.5: without censoring (April)  #################
###############################################################

setwd("/Users/Redelmeier/Dropbox/Research/Simulations/Feb 6 New Code")

source("/Users/Redelmeier/Dropbox/Research/Simulations/Feb 6 New Code/Source_keep_max.R")

varList <- 
  varlist(
    n.sim = list(type = "N", expr = quote(N[sim]), value = 1000),
    n = list(type = "frozen",  value = 1000),
    num_inj = list(type = "frozen", value = 1000),
    xi = list(type = "grid", expr = quote(xi),   value = c(0.1, 0.5, 1.5)),
    rate_exp = list(type = "frozen",  value = 1),
    ne = list(type = "frozen",  value = 50),
    method = list(type = "inner", expr = quote(method), value = c("MLE", "CensMLE")) 
  )


doOne <- function(n, xi, rate_exp, num_inj, method, ne){
  stopifnot(require(QRM), require(dplyr), sapply(list(n, xi, rate_exp, num_inj, ne), is.numeric))
  met <- rep(NA, length(method))
  names(met) <- method
  
  full <- gen_data_nocens(xi, n, num_inj, rate_exp)
  
  m <- maximum(full=full)
  
  exc <- excess(data = m, ne = ne)
  
  for(i in 1:length(method)){
    met[i] <-  mle(exc, method = method[i])
  }
  met
}

# Parallel computing
no_cores <- detectCores()-1
cl <- makeCluster(no_cores)
clusterExport(cl,"df")
clusterEvalQ(cl,library(QRM))
clusterEvalQ(cl,source("/Users/Redelmeier/Dropbox/Research/Simulations/Feb 6 New Code/Source_keep_max.R"))

#res <- doForeach(varList, cluster=cl, seed="seq", doOne=doOne,  sfile="max_no_cens60.rds", monitor=TRUE)
res <- doForeach(varList, cluster=cl, seed="seq", doOne=doOne,  sfile="max_no_cens100.rds", monitor=TRUE)
val <- getArray(res)
df <- array2df(val)

df_MLE_1 <- df %>% dplyr::filter(method == "MLE") %>% dplyr::filter(xi == 0.1) %>% dplyr::summarise(bias = mean(value) - 0.1, sd = sd(value)) %>% 
  dplyr::mutate(MSE = bias^2 + sd^2)
df_MLE_5 <- df %>% dplyr::filter(method == "MLE") %>% dplyr::filter(xi == 0.5) %>% dplyr::summarise(bias = mean(value) - 0.5, sd = sd(value))%>% 
  dplyr::mutate(MSE = bias^2 + sd^2)
df_MLE_15 <- df %>% dplyr::filter(method == "MLE") %>% dplyr::filter(xi == 1.5) %>% dplyr::summarise(bias = mean(value) - 1.5, sd = sd(value))%>% 
  dplyr::mutate(MSE = bias^2 + sd^2)
df_Cens_1 <- df %>% dplyr::filter(method == "CensMLE") %>% dplyr::filter(xi == 0.1) %>% dplyr::summarise(bias = mean(value) - 0.1, sd = sd(value))%>% 
  dplyr::mutate(MSE = bias^2 + sd^2)
df_Cens_5 <- df %>% dplyr::filter(method == "CensMLE") %>% dplyr::filter(xi == 0.5) %>% dplyr::summarise(bias = mean(value) - 0.5, sd = sd(value))%>% 
  dplyr::mutate(MSE = bias^2 + sd^2)
df_Cens_15 <- df %>% dplyr::filter(method == "CensMLE") %>% dplyr::filter(xi == 1.5) %>% dplyr::summarise(bias = mean(value) - 1.5, sd = sd(value))%>% 
  dplyr::mutate(MSE = bias^2 + sd^2)

MSE_NoCens <- matrix(NA, nrow = 3, ncol = 4)
colnames(MSE_NoCens) <- c("xi", "Bias", "sd", "MSE")
MSE_NoCens[,1] <- c(0.1, 0.5, 1.5)
MSE_NoCens[1,2] <- df_MLE_1[[1]]
MSE_NoCens[2,2] <- df_MLE_5[[1]]
MSE_NoCens[3,2] <- df_MLE_15[[1]]
MSE_NoCens[1,3] <- df_MLE_1[[2]]
MSE_NoCens[2,3] <- df_MLE_5[[2]]
MSE_NoCens[3,3] <- df_MLE_15[[2]]
MSE_NoCens[1,4] <- df_MLE_1[[3]]
MSE_NoCens[2,4] <- df_MLE_5[[3]]
MSE_NoCens[3,4] <- df_MLE_15[[3]]


xtable(MSE_NoCens, digits = c(1, 1, 3, 3, 3))


#setwd("/Users/Redelmeier/Dropbox/Research/Presentations/To_Show_Plots")

setwd("/Users/Redelmeier/Dropbox/Research/Thesis")

pdf("Sim_March_No_Cens.pdf",width=5, height=5)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 5, 0, 0))
boxplot(I(value-as.numeric(as.character(df$xi)))~xi, data = df)
abline(h=0,col="red")
mtext(text="True xi", side=1, line=4, las=1)
mtext(text="Observed Bias", side=2, line=4, las=0)
dev.off()


###############################################################
#### Study 1.5: without censoring (April - change xi)  ########
###############################################################

setwd("/Users/Redelmeier/Dropbox/Research/Simulations/Feb 6 New Code")

source("/Users/Redelmeier/Dropbox/Research/Simulations/Feb 6 New Code/Source_keep_max.R")

varList <- 
  varlist(
    n.sim = list(type = "N", expr = quote(N[sim]), value = 3000),
    n = list(type = "frozen",  value = 1000),
    num_inj = list(type = "frozen", value = 60),
    xi = list(type = "grid", expr = quote(xi),   value = c(0.1, 0.3, 0.5)),
    rate_exp = list(type = "frozen",  value = 1),
    ne = list(type = "frozen",  value = 150),
    method = list(type = "inner", expr = quote(method), value = c("MLE", "CensMLE")) 
  )


doOne <- function(n, xi, rate_exp, num_inj, method, ne){
  stopifnot(require(QRM), require(dplyr), sapply(list(n, xi, rate_exp, num_inj, ne), is.numeric))
  met <- rep(NA, length(method))
  names(met) <- method
  
  full <- gen_data_nocens(xi, n, num_inj, rate_exp)
  
  m <- maximum(full=full)
  
  exc <- excess(data = m, ne = ne)
  
  for(i in 1:length(method)){
    met[i] <-  mle(exc, method = method[i])
  }
  met
}

# Parallel computing
no_cores <- detectCores()-1
cl <- makeCluster(no_cores)
clusterExport(cl,"df")
clusterEvalQ(cl,library(QRM))
clusterEvalQ(cl,source("/Users/Redelmeier/Dropbox/Research/Simulations/Feb 6 New Code/Source_keep_max.R"))

#res <- doForeach(varList, cluster=cl, seed="seq", doOne=doOne,  sfile="max_no_cens60_05_3000.rds", monitor=TRUE)
res <- doForeach(varList, cluster=cl, seed="seq", doOne=doOne,  sfile="max_no_cens60_05.rds", monitor=TRUE)
val <- getArray(res)
df <- array2df(val)

df_MLE_1 <- df %>% dplyr::filter(method == "MLE") %>% dplyr::filter(xi == 0.1) %>% dplyr::summarise(bias = mean(value) - 0.1, sd = sd(value)) %>% 
  dplyr::mutate(MSE = bias^2 + sd^2)
df_MLE_5 <- df %>% dplyr::filter(method == "MLE") %>% dplyr::filter(xi == 0.3) %>% dplyr::summarise(bias = mean(value) - 0.3, sd = sd(value))%>% 
  dplyr::mutate(MSE = bias^2 + sd^2)
df_MLE_15 <- df %>% dplyr::filter(method == "MLE") %>% dplyr::filter(xi == 0.5) %>% dplyr::summarise(bias = mean(value) - 0.5, sd = sd(value))%>% 
  dplyr::mutate(MSE = bias^2 + sd^2)
df_Cens_1 <- df %>% dplyr::filter(method == "CensMLE") %>% dplyr::filter(xi == 0.1) %>% dplyr::summarise(bias = mean(value) - 0.1, sd = sd(value))%>% 
  dplyr::mutate(MSE = bias^2 + sd^2)
df_Cens_5 <- df %>% dplyr::filter(method == "CensMLE") %>% dplyr::filter(xi == 0.3) %>% dplyr::summarise(bias = mean(value) - 0.3, sd = sd(value))%>% 
  dplyr::mutate(MSE = bias^2 + sd^2)
df_Cens_15 <- df %>% dplyr::filter(method == "CensMLE") %>% dplyr::filter(xi == 0.5) %>% dplyr::summarise(bias = mean(value) - 0.5, sd = sd(value))%>% 
  dplyr::mutate(MSE = bias^2 + sd^2)

MSE_NoCens <- matrix(NA, nrow = 3, ncol = 4)
colnames(MSE_NoCens) <- c("xi", "Bias", "sd", "MSE")
MSE_NoCens[,1] <- c(0.1, 0.3, 0.5)
MSE_NoCens[1,2] <- df_MLE_1[[1]]
MSE_NoCens[2,2] <- df_MLE_5[[1]]
MSE_NoCens[3,2] <- df_MLE_15[[1]]
MSE_NoCens[1,3] <- df_MLE_1[[2]]
MSE_NoCens[2,3] <- df_MLE_5[[2]]
MSE_NoCens[3,3] <- df_MLE_15[[2]]
MSE_NoCens[1,4] <- df_MLE_1[[3]]
MSE_NoCens[2,4] <- df_MLE_5[[3]]
MSE_NoCens[3,4] <- df_MLE_15[[3]]


xtable(MSE_NoCens, digits = c(1, 1, 3, 3, 3))


#setwd("/Users/Redelmeier/Dropbox/Research/Presentations/To_Show_Plots")

setwd("/Users/Redelmeier/Dropbox/Research/Thesis")

pdf("Sim_March_No_Cens.pdf",width=5, height=5)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 5, 0, 0))
boxplot(I(value-as.numeric(as.character(df$xi)))~xi, data = df)
abline(h=0,col="red")
mtext(text="True xi", side=1, line=4, las=1)
mtext(text="Observed Bias", side=2, line=4, las=0)
dev.off()



###############################################################
####### Study 2: with censoring  ##############################
###############################################################

setwd("/Users/Redelmeier/Dropbox/Research/Simulations/Feb 6 New Code")

source("/Users/Redelmeier/Dropbox/Research/Simulations/Feb 6 New Code/Source_keep_max.R")

varList2 <- 
  varlist(
    n.sim = list(type = "N", expr = quote(N[sim]), value = 1000),
    n = list(type = "frozen",  value = 1000),
    num_inj = list(type = "frozen", value = 60),
    xi = list(type = "frozen", value = 0.5),
    censor = list(type = "grid", expr = quote(censor), value = c(36, 60, 90)),
    rate_exp = list(type = "frozen",  value = 1),
    ne = list(type = "frozen",  value = 50),
    method = list(type = "inner", expr = quote(method), value = c("MLE", "CensMLE")) 
  )


doOne2 <- function(n, xi, rate_exp, censor, num_inj, method, ne){
  stopifnot(require(QRM), require(dplyr), sapply(list(n, xi, censor, rate_exp, num_inj, ne), is.numeric))
  met <- rep(NA, length(method))
  names(met) <- method
  
  full <- gen_data(xi=xi, censor = censor, n=n, num_inj=num_inj, rate_exp=rate_exp)
  
  m <- maximum(full=full)
  
  exc <- excess(data = m, ne = ne)
  
  for(i in 1:length(method)){
    met[i] <-  mle(exc, method = method[i])
  }
  met
}

# Parallel computing
no_cores <- detectCores()-1
cl <- makeCluster(no_cores)
clusterExport(cl,"df")
clusterEvalQ(cl,library(QRM))
clusterEvalQ(cl,source("/Users/Redelmeier/Dropbox/Research/Simulations/Feb 6 New Code/Source_keep_max.R"))

res2 <- doForeach(varList2, cluster=cl, seed="seq", doOne=doOne2,  sfile="max_cens.rds", monitor=TRUE)
val2 <- getArray(res2)
df2 <- array2df(val2)


df2_MLE_36 <- df2 %>% dplyr::filter(method == "MLE") %>% dplyr::filter(censor == 36) %>% dplyr::summarise(bias = mean(value) - 0.5, sd = sd(value))%>% 
  dplyr::mutate(MSE = bias^2 + sd^2)
df2_MLE_60 <- df2 %>% dplyr::filter(method == "MLE") %>% dplyr::filter(censor == 60) %>% dplyr::summarise(bias = mean(value) - 0.5, sd = sd(value))%>% 
  dplyr::mutate(MSE = bias^2 + sd^2)
df2_MLE_90 <- df2 %>% dplyr::filter(method == "MLE") %>% dplyr::filter(censor == 90) %>% dplyr::summarise(bias = mean(value) - 0.5, sd = sd(value))%>% 
  dplyr::mutate(MSE = bias^2 + sd^2)
df2_Cens_36 <- df2 %>% dplyr::filter(method == "CensMLE") %>% dplyr::filter(censor == 36) %>% dplyr::summarise(bias = mean(value) - 0.5, sd = sd(value))%>% 
  dplyr::mutate(MSE = bias^2 + sd^2)
df2_Cens_60 <- df2 %>% dplyr::filter(method == "CensMLE") %>% dplyr::filter(censor == 60) %>% dplyr::summarise(bias = mean(value) - 0.5, sd = sd(value))%>% 
  dplyr::mutate(MSE = bias^2 + sd^2)
df2_Cens_90 <- df2 %>% dplyr::filter(method == "CensMLE") %>% dplyr::filter(censor == 90) %>% dplyr::summarise(bias = mean(value) - 0.5, sd = sd(value))%>% 
  dplyr::mutate(MSE = bias^2 + sd^2)



MSE_Cens <- matrix(NA, nrow = 3, ncol = 7)
colnames(MSE_Cens) <- c("C", "Bias", "sd", "MSE", "Bias", "sd", "MSE")
MSE_Cens[,1] <- c(36, 60, 90)
MSE_Cens[1,2] <- df2_MLE_36[[1]]
MSE_Cens[1,3] <- df2_MLE_36[[2]]
MSE_Cens[1,4] <- df2_MLE_36[[3]]

MSE_Cens[2,2] <- df2_MLE_60[[1]]
MSE_Cens[2,3] <- df2_MLE_60[[2]]
MSE_Cens[2,4] <- df2_MLE_60[[3]]

MSE_Cens[3,2] <- df2_MLE_90[[1]]
MSE_Cens[3,3] <- df2_MLE_90[[2]]
MSE_Cens[3,4] <- df2_MLE_90[[3]]

MSE_Cens[1,5] <- df2_Cens_36[[1]]
MSE_Cens[1,6] <- df2_Cens_36[[2]]
MSE_Cens[1,7] <- df2_Cens_36[[3]]

MSE_Cens[2,5] <- df2_Cens_60[[1]]
MSE_Cens[2,6] <- df2_Cens_60[[2]]
MSE_Cens[2,7] <- df2_Cens_60[[3]]

MSE_Cens[3,5] <- df2_Cens_90[[1]]
MSE_Cens[3,6] <- df2_Cens_90[[2]]
MSE_Cens[3,7] <- df2_Cens_90[[3]]


xtable(MSE_Cens, digits = c(1, 1, 3, 3, 3, 3, 3, 3))


#setwd("/Users/Redelmeier/Dropbox/Research/Presentations/To_Show_Plots")

setwd("/Users/Redelmeier/Dropbox/Research/Thesis")

pdf("Sim_March.pdf",width=6, height=5)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 5, 0, 0))
boxplot(I(value-0.5)~ censor+method, data = df2, xaxt = "n", main= paste("Non-Censored MLE"," "," Censored MLE "))
abline(h=0,col="red")
abline(v=3.5,col="grey")
axis(side=1, at= 1:6, labels= c("36", "60", "90", "36", "60", "90"))
mtext(text="End-of-Followup Time", side=1, line=4, las=1)
mtext(text="Observed Bias", side=2, line=4, las=0)
dev.off()





###############################################################
####### Study 3: count censored  ##############################
###############################################################

setwd("/Users/Redelmeier/Dropbox/Research/Simulations/Feb 6 New Code")

source("/Users/Redelmeier/Dropbox/Research/Simulations/Feb 6 New Code/Source_keep_max.R")

varList3 <- 
  varlist(
    n.sim = list(type = "N", expr = quote(N[sim]), value = 1000),
    n = list(type = "frozen",  value = 1000),
    num_inj = list(type = "frozen", value = 60),
    xi = list(type = "frozen", value = 0.5),
    censor = list(type = "grid", expr = quote(censor), value = c(36, 60, 90)),
    rate_exp = list(type = "frozen",  value = 1),
    ne = list(type = "frozen",  value = 50),
    method = list(type = "inner", expr = quote(method), value = c("N")) 
  )


doOne2 <- function(n, xi, rate_exp, censor, num_inj, method, ne){
  stopifnot(require(QRM), require(dplyr), sapply(list(n, xi, censor, rate_exp, num_inj, ne), is.numeric))
  met <- rep(NA, length(method))
  names(met) <- method
  
  full <- gen_data(xi=xi, censor = censor, n=n, num_inj=num_inj, rate_exp=rate_exp)
  
  m <- maximum(full=full)
  
  exc <- excess(data = m, ne = ne)
  
  for(i in 1:length(method)){
    met[i] <-  mle(exc, method = method[i])
  }
  met
}

# Parallel computing
no_cores <- detectCores()-1
cl <- makeCluster(no_cores)
clusterExport(cl,"df")
clusterEvalQ(cl,library(QRM))
clusterEvalQ(cl,source("/Users/Redelmeier/Dropbox/Research/Simulations/Feb 6 New Code/Source_keep_max.R"))

res3 <- doForeach(varList3, cluster=cl, seed="seq", doOne=doOne2,  sfile="max_cens_N.rds", monitor=TRUE)
val3 <- getArray(res3)
df3 <- array2df(val3)

#setwd("/Users/Redelmeier/Dropbox/Research/Presentations/To_Show_Plots")

setwd("/Users/Redelmeier/Dropbox/Research/Thesis")

pdf("Sim_March_N.pdf",width=5, height=5)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 5, 0, 0))
boxplot(I(value) ~ censor, data = df3, xaxt = "n", yaxt = "n")
axis(side=1, at= 1:3, labels= c("36", "60", "120"))
axis(side=2, at= c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7), labels= c("20", "30", "40", "50", "60", "70" ))
mtext(text="End-of-Followup Time", side=1, line=4, las=1)
mtext(text="Percentage of Censored Observations", side=2, line=4, las=0)
dev.off()





