setwd("/Users/Redelmeier/Dropbox/Research/Simulations/Feb 6 New Code")

library(simsalapar)
library(parallel)
library(foreach)
library(doParallel)
library(QRM)
library(dplyr)
library(xtable)

###############################################################
######### Study 1: Delete Censored Obs  #######################
###############################################################


source("/Users/Redelmeier/Dropbox/Research/Simulations/Feb 6 New Code/Source_delete.R")

varList <- 
  varlist(
    n.sim = list(type = "N", expr = quote(N[sim]), value = 1000),
    n = list(type = "frozen", value = 1000),
    num_inj = list(type = "frozen", value = 30),
    xi = list(type = "frozen", value = 0.5),
    censor = list(type = "grid", expr = quote(censor), value = c(3, 6, 12, 36)),
    beta = list(type = "frozen",  value = 1),
    rate_exp = list(type = "frozen", value = 1),
    method = list(type = "inner", expr = quote(method), value = c("MLE", "CensMLE")) 
  )

doOne <- function(n, xi, censor, beta, rate_exp, num_inj, method){
  stopifnot(require(QRM), require(dplyr), sapply(list(n, xi, censor, rate_exp, num_inj, beta), is.numeric))
  met <- rep(NA, length(method))
  names(met) <- method
  
  full <- gen_data(censor, xi, n, num_inj, beta, rate_exp)
  
  for(i in 1:length(method)){
    met[i] <-  mle(full, method = method[i])
  }
  met
}

# Parallel computing
no_cores <- detectCores()-1
cl <- makeCluster(no_cores)
clusterExport(cl,"df")
clusterEvalQ(cl,library(QRM))
clusterEvalQ(cl,source("/Users/Redelmeier/Dropbox/Research/Simulations/Feb 6 New Code/Source_delete.R"))

res <- doForeach(varList,cluster=cl, seed="seq", doOne=doOne,  sfile="Cens_remove.rds",monitor=TRUE)
val <- getArray(res)
df <- array2df(val)

setwd("/Users/Redelmeier/Dropbox/Research/Presentations/To_Show_Plots")

pdf("Sim2_delete_March.pdf",width=5, height=5)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 5, 0, 0))
boxplot(I(value-0.5)~censor + method, data = df, xaxt = "n", main = "MLE/CensMLE")
abline(h=0,col="red")
abline(v=4.5,col="grey")
axis(side=1, at= 1:8, label = c("3", "6", "12", "36","3", "6", "12", "36"))
mtext(text="Censoring Time", side=1, line=4, las=1)
mtext(text="Observed Bias", side=2, line=4, las=0)
dev.off()


###############################################################
######### Study 2: Keep Censored Obs  #########################
###############################################################

setwd("/Users/Redelmeier/Dropbox/Research/Simulations/Feb 6 New Code")

source("/Users/Redelmeier/Dropbox/Research/Simulations/Feb 6 New Code/Source_keep.R")

varList <- 
  varlist(
    n.sim = list(type = "N", expr = quote(N[sim]), value = 1000),
    n = list(type = "frozen", value = 1000),
    num_inj = list(type = "frozen", value = 30),
    xi = list(type = "frozen", value = 0.5),
    censor = list(type = "grid", expr = quote(censor), value = c(3, 6, 12, 36)),
    beta = list(type = "frozen",  value = 1),
    rate_exp = list(type = "frozen", value = 1),
    method = list(type = "inner", expr = quote(method), value = c("MLE", "CensMLE")) 
  )

doOne <- function(n, xi, censor, beta, rate_exp, num_inj, method){
  stopifnot(require(QRM), require(dplyr), sapply(list(n, xi, censor, rate_exp, num_inj, beta), is.numeric))
  met <- rep(NA, length(method))
  names(met) <- method
  
  full <- gen_data(censor, xi, n, num_inj, beta, rate_exp)
  
  for(i in 1:length(method)){
    met[i] <-  mle(full, method = method[i])
  }
  met
}

# Parallel computing
no_cores <- detectCores()-1
cl <- makeCluster(no_cores)
clusterExport(cl,"df")
clusterEvalQ(cl,library(QRM))
clusterEvalQ(cl,source("/Users/Redelmeier/Dropbox/Research/Simulations/Feb 6 New Code/Source_keep.R"))

res2 <- doForeach(varList,cluster=cl, seed="seq", doOne=doOne,  sfile="Cens_keep.rds",monitor=TRUE)
val2 <- getArray(res2)
df2 <- array2df(val2)

setwd("/Users/Redelmeier/Dropbox/Research/Presentations/To_Show_Plots")

pdf("Sim2_keep_March.pdf",width=5, height=5)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 5, 0, 0))
boxplot(I(value-0.5)~censor+method, data = df2, xaxt = "n")
abline(h=0,col="red")
abline(v=4.5,col="grey")
axis(side=1, at= 1:8, labels= c("3", "6", "12", "36", "3", "6", "12", "36"))
mtext(text="Censoring Time", side=1, line=4, las=0)
mtext(text="Observed Bias", side=2, line=4, las=0)
dev.off()

###############################################################
####### Study 2 Extra: Both Keep and Delete (all data)  #######
###############################################################

setwd("/Users/Redelmeier/Dropbox/Research/Simulations/Feb 6 New Code")

source("/Users/Redelmeier/Dropbox/Research/Simulations/Feb 6 New Code/Source_both.R")

varList_all <- 
  varlist(
    n.sim = list(type = "N", expr = quote(N[sim]), value = 1000),
    n = list(type = "frozen", value = 1000),
    num_inj = list(type = "frozen", value = 30),
    xi = list(type = "frozen", value = 0.5),
    censor = list(type = "grid", expr = quote(censor), value = c(3, 6, 12, 36)),
    beta = list(type = "frozen",  value = 1),
    rate_exp = list(type = "frozen", value = 1),
    method = list(type = "inner", expr = quote(method), value = c("MLE_All", "MLE_Delete", "CensMLE")) 
  )

doOne <- function(n, xi, censor, beta, rate_exp, num_inj, method){
  stopifnot(require(QRM), require(dplyr), sapply(list(n, xi, censor, rate_exp, num_inj, beta), is.numeric))
  met <- rep(NA, length(method))
  names(met) <- method
  
  full <- gen_data(censor, xi, n, num_inj, beta, rate_exp)
  
  for(i in 1:length(method)){
    met[i] <-  mle(full, method = method[i])
  }
  met
}

# Parallel computing
no_cores <- detectCores()-1
cl <- makeCluster(no_cores)
clusterExport(cl,"df")
clusterEvalQ(cl,library(QRM))
clusterEvalQ(cl,source("/Users/Redelmeier/Dropbox/Research/Simulations/Feb 6 New Code/Source_both.R"))

res_all <- doForeach(varList_all, cluster=cl, seed="seq", doOne=doOne,  sfile="Cens_keep_all.rds", monitor=TRUE)
val_all <- getArray(res_all)
df_all <- array2df(val_all)

#setwd("/Users/Redelmeier/Dropbox/Research/Presentations/To_Show_Plots")
setwd("/Users/Redelmeier/Dropbox/Research/Thesis")

pdf("Sim2_all.pdf",width=7, height=5)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 5, 0, 0))
boxplot(I(value-0.5)~censor+method, data = df_all, xaxt = "n",main= paste("Non-Censored MLE "," "," Delete Obs"," "," Censored MLE"))
abline(h=0,col="red")
abline(v=4.5,col="grey")
abline(v=8.5,col="grey")
axis(side=1, at= 1:12, labels= c("3", "6", "12", "36","3", "6", "12", "36", "3", "6", "12", "36"))
mtext(text="End-Of-Followup Time", side=1, line=4, las=0)
mtext(text="Observed Bias", side=2, line=4, las=0)
dev.off()




df_all_MLE_3 <- df_all %>% dplyr::filter(method == "MLE_All") %>% dplyr::filter(censor == 3) %>% dplyr::summarise(bias = mean(value) - 0.5, sd = sd(value))%>% 
  dplyr::mutate(MSE = bias^2 + sd^2)
df_all_MLE_6 <- df_all %>% dplyr::filter(method == "MLE_All") %>% dplyr::filter(censor == 6) %>% dplyr::summarise(bias = mean(value) - 0.5, sd = sd(value))%>% 
  dplyr::mutate(MSE = bias^2 + sd^2)
df_all_MLE_12 <- df_all %>% dplyr::filter(method == "MLE_All") %>% dplyr::filter(censor == 12) %>% dplyr::summarise(bias = mean(value) - 0.5, sd = sd(value))%>% 
  dplyr::mutate(MSE = bias^2 + sd^2)
df_all_MLE_36 <- df_all %>% dplyr::filter(method == "MLE_All") %>% dplyr::filter(censor == 36) %>% dplyr::summarise(bias = mean(value) - 0.5, sd = sd(value))%>% 
  dplyr::mutate(MSE = bias^2 + sd^2)

df_all_DelMLE_3 <- df_all %>% dplyr::filter(method == "MLE_Delete") %>% dplyr::filter(censor == 3) %>% dplyr::summarise(bias = mean(value) - 0.5, sd = sd(value))%>% 
  dplyr::mutate(MSE = bias^2 + sd^2)
df_all_DelMLE_6 <- df_all %>% dplyr::filter(method == "MLE_Delete") %>% dplyr::filter(censor == 6) %>% dplyr::summarise(bias = mean(value) - 0.5, sd = sd(value))%>% 
  dplyr::mutate(MSE = bias^2 + sd^2)
df_all_DelMLE_12 <- df_all %>% dplyr::filter(method == "MLE_Delete") %>% dplyr::filter(censor == 12) %>% dplyr::summarise(bias = mean(value) - 0.5, sd = sd(value))%>% 
  dplyr::mutate(MSE = bias^2 + sd^2)
df_all_DelMLE_36 <- df_all %>% dplyr::filter(method == "MLE_Delete") %>% dplyr::filter(censor == 36) %>% dplyr::summarise(bias = mean(value) - 0.5, sd = sd(value))%>% 
  dplyr::mutate(MSE = bias^2 + sd^2)

df_all_CensMLE_3 <- df_all %>% dplyr::filter(method == "CensMLE") %>% dplyr::filter(censor == 3) %>% dplyr::summarise(bias = mean(value) - 0.5, sd = sd(value))%>% 
  dplyr::mutate(MSE = bias^2 + sd^2)
df_all_CensMLE_6 <- df_all %>% dplyr::filter(method == "CensMLE") %>% dplyr::filter(censor == 6) %>% dplyr::summarise(bias = mean(value) - 0.5, sd = sd(value))%>% 
  dplyr::mutate(MSE = bias^2 + sd^2)
df_all_CensMLE_12 <- df_all %>% dplyr::filter(method == "CensMLE") %>% dplyr::filter(censor == 12) %>% dplyr::summarise(bias = mean(value) - 0.5, sd = sd(value))%>% 
  dplyr::mutate(MSE = bias^2 + sd^2)
df_all_CensMLE_36 <- df_all %>% dplyr::filter(method == "CensMLE") %>% dplyr::filter(censor == 36) %>% dplyr::summarise(bias = mean(value) - 0.5, sd = sd(value))%>% 
  dplyr::mutate(MSE = bias^2 + sd^2)


MSE_Cens <- matrix(NA, nrow = 4, ncol = 9)
colnames(MSE_Cens) <- c("Bias", "sd", "MSE", "Bias", "sd", "MSE", "Bias", "sd", "MSE")
rownames(MSE_Cens) <- c(3, 6, 12, 26)

MSE_Cens[1,1] <- df_all_MLE_3[[1]]
MSE_Cens[1,2] <- df_all_MLE_3[[2]]
MSE_Cens[1,3] <- df_all_MLE_3[[3]]

MSE_Cens[2,1] <- df_all_MLE_6[[1]]
MSE_Cens[2,2] <- df_all_MLE_6[[2]]
MSE_Cens[2,3] <- df_all_MLE_6[[3]]

MSE_Cens[3,1] <- df_all_MLE_12[[1]]
MSE_Cens[3,2] <- df_all_MLE_12[[2]]
MSE_Cens[3,3] <- df_all_MLE_12[[3]]

MSE_Cens[4,1] <- df_all_MLE_36[[1]]
MSE_Cens[4,2] <- df_all_MLE_36[[2]]
MSE_Cens[4,3] <- df_all_MLE_36[[3]]

MSE_Cens[1,4] <- df_all_DelMLE_3[[1]]
MSE_Cens[1,5] <- df_all_DelMLE_3[[2]]
MSE_Cens[1,6] <- df_all_DelMLE_3[[3]]

MSE_Cens[2,4] <- df_all_DelMLE_6[[1]]
MSE_Cens[2,5] <- df_all_DelMLE_6[[2]]
MSE_Cens[2,6] <- df_all_DelMLE_6[[3]]

MSE_Cens[3,4] <- df_all_DelMLE_12[[1]]
MSE_Cens[3,5] <- df_all_DelMLE_12[[2]]
MSE_Cens[3,6] <- df_all_DelMLE_12[[3]]

MSE_Cens[4,4] <- df_all_DelMLE_36[[1]]
MSE_Cens[4,5] <- df_all_DelMLE_36[[2]]
MSE_Cens[4,6] <- df_all_DelMLE_36[[3]]

MSE_Cens[1,7] <- df_all_CensMLE_3[[1]]
MSE_Cens[1,8] <- df_all_CensMLE_3[[2]]
MSE_Cens[1,9] <- df_all_CensMLE_3[[3]]

MSE_Cens[2,7] <- df_all_CensMLE_6[[1]]
MSE_Cens[2,8] <- df_all_CensMLE_6[[2]]
MSE_Cens[2,9] <- df_all_CensMLE_6[[3]]

MSE_Cens[3,7] <- df_all_CensMLE_12[[1]]
MSE_Cens[3,8] <- df_all_CensMLE_12[[2]]
MSE_Cens[3,9] <- df_all_CensMLE_12[[3]]

MSE_Cens[4,7] <- df_all_CensMLE_36[[1]]
MSE_Cens[4,8] <- df_all_CensMLE_36[[2]]
MSE_Cens[4,9] <- df_all_CensMLE_36[[3]]



xtable(MSE_Cens, digits = c(3, 3, 3, 3, 3, 3, 3,3,3,3))

###############################################################
####### Study 2.5 Extra: Both Keep and Delete (all data)  #######
###############################################################

setwd("/Users/Redelmeier/Dropbox/Research/Simulations/Feb 6 New Code")

source("/Users/Redelmeier/Dropbox/Research/Simulations/Feb 6 New Code/Source_both.R")

varList_all <- 
  varlist(
    n.sim = list(type = "N", expr = quote(N[sim]), value = 1000),
    n = list(type = "frozen", value = 1000),
    num_inj = list(type = "frozen", value = 60),
    xi = list(type = "frozen", value = 0.5),
    censor = list(type = "grid", expr = quote(censor), value = c(36, 60, 90)),
    beta = list(type = "frozen",  value = 1),
    rate_exp = list(type = "frozen", value = 1),
    method = list(type = "inner", expr = quote(method), value = c("MLE_All", "MLE_Delete", "CensMLE")) 
  )

doOne <- function(n, xi, censor, beta, rate_exp, num_inj, method){
  stopifnot(require(QRM), require(dplyr), sapply(list(n, xi, censor, rate_exp, num_inj, beta), is.numeric))
  met <- rep(NA, length(method))
  names(met) <- method
  
  full <- gen_data(censor, xi, n, num_inj, beta, rate_exp)
  
  for(i in 1:length(method)){
    met[i] <-  mle(full, method = method[i])
  }
  met
}

# Parallel computing
no_cores <- detectCores()-1
cl <- makeCluster(no_cores)
clusterExport(cl,"df")
clusterEvalQ(cl,library(QRM))
clusterEvalQ(cl,source("/Users/Redelmeier/Dropbox/Research/Simulations/Feb 6 New Code/Source_both.R"))

res_all <- doForeach(varList_all, cluster=cl, seed="seq", doOne=doOne,  sfile="Cens_keep_all60.rds", monitor=TRUE)
val_all <- getArray(res_all)
df_all <- array2df(val_all)

#setwd("/Users/Redelmeier/Dropbox/Research/Presentations/To_Show_Plots")
setwd("/Users/Redelmeier/Dropbox/Research/Thesis")

pdf("Sim2_all.pdf",width=7, height=5)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 5, 0, 0))
boxplot(I(value-0.5)~censor+method, data = df_all, xaxt = "n",main= paste("Non-Censored MLE   "," ","    Delete Obs "," ","   Censored MLE"))
abline(h=0,col="red")
abline(v=3.5,col="grey")
abline(v=6.5,col="grey")
axis(side=1, at= 1:12, labels= c("36", "60", "90", "36", "60", "90", "36", "60", "90"))
mtext(text="End-Of-Followup Time", side=1, line=4, las=0)
mtext(text="Observed Bias", side=2, line=4, las=0)
dev.off()




df_all_MLE_3 <- df_all %>% dplyr::filter(method == "MLE_All") %>% dplyr::filter(censor == 36) %>% dplyr::summarise(bias = mean(value) - 0.5, sd = sd(value))%>% 
  dplyr::mutate(MSE = bias^2 + sd^2)
df_all_MLE_6 <- df_all %>% dplyr::filter(method == "MLE_All") %>% dplyr::filter(censor == 60) %>% dplyr::summarise(bias = mean(value) - 0.5, sd = sd(value))%>% 
  dplyr::mutate(MSE = bias^2 + sd^2)
df_all_MLE_12 <- df_all %>% dplyr::filter(method == "MLE_All") %>% dplyr::filter(censor == 90) %>% dplyr::summarise(bias = mean(value) - 0.5, sd = sd(value))%>% 
  dplyr::mutate(MSE = bias^2 + sd^2)


df_all_DelMLE_3 <- df_all %>% dplyr::filter(method == "MLE_Delete") %>% dplyr::filter(censor == 36) %>% dplyr::summarise(bias = mean(value) - 0.5, sd = sd(value))%>% 
  dplyr::mutate(MSE = bias^2 + sd^2)
df_all_DelMLE_6 <- df_all %>% dplyr::filter(method == "MLE_Delete") %>% dplyr::filter(censor == 60) %>% dplyr::summarise(bias = mean(value) - 0.5, sd = sd(value))%>% 
  dplyr::mutate(MSE = bias^2 + sd^2)
df_all_DelMLE_12 <- df_all %>% dplyr::filter(method == "MLE_Delete") %>% dplyr::filter(censor == 90) %>% dplyr::summarise(bias = mean(value) - 0.5, sd = sd(value))%>% 
  dplyr::mutate(MSE = bias^2 + sd^2)


df_all_CensMLE_3 <- df_all %>% dplyr::filter(method == "CensMLE") %>% dplyr::filter(censor == 36) %>% dplyr::summarise(bias = mean(value) - 0.5, sd = sd(value))%>% 
  dplyr::mutate(MSE = bias^2 + sd^2)
df_all_CensMLE_6 <- df_all %>% dplyr::filter(method == "CensMLE") %>% dplyr::filter(censor == 60) %>% dplyr::summarise(bias = mean(value) - 0.5, sd = sd(value))%>% 
  dplyr::mutate(MSE = bias^2 + sd^2)
df_all_CensMLE_12 <- df_all %>% dplyr::filter(method == "CensMLE") %>% dplyr::filter(censor == 90) %>% dplyr::summarise(bias = mean(value) - 0.5, sd = sd(value))%>% 
  dplyr::mutate(MSE = bias^2 + sd^2)



MSE_Cens <- matrix(NA, nrow = 3, ncol = 9)
colnames(MSE_Cens) <- c("Bias", "sd", "MSE", "Bias", "sd", "MSE", "Bias", "sd", "MSE")
rownames(MSE_Cens) <- c(36, 60, 90)

MSE_Cens[1,1] <- df_all_MLE_3[[1]]
MSE_Cens[1,2] <- df_all_MLE_3[[2]]
MSE_Cens[1,3] <- df_all_MLE_3[[3]]

MSE_Cens[2,1] <- df_all_MLE_6[[1]]
MSE_Cens[2,2] <- df_all_MLE_6[[2]]
MSE_Cens[2,3] <- df_all_MLE_6[[3]]

MSE_Cens[3,1] <- df_all_MLE_12[[1]]
MSE_Cens[3,2] <- df_all_MLE_12[[2]]
MSE_Cens[3,3] <- df_all_MLE_12[[3]]


MSE_Cens[1,4] <- df_all_DelMLE_3[[1]]
MSE_Cens[1,5] <- df_all_DelMLE_3[[2]]
MSE_Cens[1,6] <- df_all_DelMLE_3[[3]]

MSE_Cens[2,4] <- df_all_DelMLE_6[[1]]
MSE_Cens[2,5] <- df_all_DelMLE_6[[2]]
MSE_Cens[2,6] <- df_all_DelMLE_6[[3]]

MSE_Cens[3,4] <- df_all_DelMLE_12[[1]]
MSE_Cens[3,5] <- df_all_DelMLE_12[[2]]
MSE_Cens[3,6] <- df_all_DelMLE_12[[3]]


MSE_Cens[1,7] <- df_all_CensMLE_3[[1]]
MSE_Cens[1,8] <- df_all_CensMLE_3[[2]]
MSE_Cens[1,9] <- df_all_CensMLE_3[[3]]

MSE_Cens[2,7] <- df_all_CensMLE_6[[1]]
MSE_Cens[2,8] <- df_all_CensMLE_6[[2]]
MSE_Cens[2,9] <- df_all_CensMLE_6[[3]]

MSE_Cens[3,7] <- df_all_CensMLE_12[[1]]
MSE_Cens[3,8] <- df_all_CensMLE_12[[2]]
MSE_Cens[3,9] <- df_all_CensMLE_12[[3]]




xtable(MSE_Cens, digits = c(3, 3, 3, 3, 3, 3, 3,3,3,3))




###############################################################
####### Study 3: Delete Censored Obs (Only max)  ##############
###############################################################

setwd("/Users/Redelmeier/Dropbox/Research/Simulations/Feb 6 New Code")

source("/Users/Redelmeier/Dropbox/Research/Simulations/Feb 6 New Code/Source_delete_max.R")

varList2 <- 
  varlist(
    n.sim = list(type = "N", expr = quote(N[sim]), value = 1000),
    n = list(type = "frozen", value = 1000),
    num_inj = list(type = "frozen", value = 30),
    xi = list(type = "grid", value = 0.5),
    censor = list(type = "grid", expr = quote(censor), value = c(3, 6, 12, 36)),
    beta = list(type = "frozen",  value = 1),
    rate_exp = list(type = "frozen", value = 1),
    ne = list(type = "frozen",  value = 50),
    method = list(type = "inner", expr = quote(method), value = c("MLE", "CensMLE")) 
  )


doOne2 <- function(n, xi, censor, beta, rate_exp, num_inj, method, ne){
  stopifnot(require(QRM), require(dplyr), sapply(list(n, xi, censor, rate_exp, num_inj, beta, ne), is.numeric))
  met <- rep(NA, length(method))
  names(met) <- method
  full <- gen_data(censor=censor, xi=xi, n=n, num_inj=num_inj, beta=beta, rate_exp=rate_exp)
  
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
clusterEvalQ(cl,source("/Users/Redelmeier/Dropbox/Research/Simulations/Feb 6 New Code/Source_delete_max.R"))


res3 <- doForeach(varList2, cluster=cl, seed="seq", doOne=doOne2,  sfile="Cens_delete_max.rds", monitor=TRUE)
val3 <- getArray(res3)
df3 <- array2df(val3)

setwd("/Users/Redelmeier/Dropbox/Research/Presentations/To_Show_Plots")

pdf("Sim2_delete_max_March.pdf",width=5, height=5)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 5, 0, 0))
boxplot(I(value-0.5)~censor + method, data = df3, xaxt = "n",main = "MLE/CensMLE")
abline(h=0,col="red")
axis(side=1, at= 1:8, labels= c("3", "6", "12","36","3", "6", "12","36"))
mtext(text="Censoring Time", side=1, line=4, las=1)
mtext(text="Observed Bias", side=2, line=4, las=0)
dev.off()


###############################################################
####### Study 4: Keep Censored Obs (Only max)  ################
###############################################################



###############################################################
####### Mini Study: Change Rate Parameter #####################
###############################################################


setwd("/Users/Redelmeier/Dropbox/Research/Simulations/Feb 6 New Code")
source("/Users/Redelmeier/Dropbox/Research/Simulations/Feb 6 New Code/Source_keep_max.R")

varListA <- 
  varlist(
    n.sim = list(type = "N", expr = quote(N[sim]), value = 1000),
    n = list(type = "frozen", value = 1000),
    num_inj = list(type = "frozen", value = 30),
    xi = list(type = "frozen", value = 0.5),
    censor = list(type = "grid", expr = quote(censor), value = c(36, 60, 90, 120)),
    beta = list(type = "frozen",  value = 1),
    rate_exp = list(type = "frozen", value = 0.2),
    ne = list(type = "frozen",  value = 50),
    method = list(type = "inner", expr = quote(method), value = c("MLE", "CensMLE")) 
  )

doOne2 <- function(n, xi, censor, beta, rate_exp, num_inj, method, ne){
  stopifnot(require(QRM), require(dplyr), sapply(list(n, xi, censor, rate_exp, num_inj, beta, ne), is.numeric))
  met <- rep(NA, length(method))
  names(met) <- method
  full <- gen_data(censor=censor, xi=xi, n=n, num_inj=num_inj, beta=beta, rate_exp=rate_exp)
  
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


resA <- doForeach(varListA, cluster=cl, seed="seq", doOne=doOne2,  sfile="Cens_keep_A.rds", monitor=TRUE)
valA <- getArray(resA)
dfA <- array2df(valA)

setwd("/Users/Redelmeier/Dropbox/Research/Presentations/To_Show_Plots")

pdf("Cens_keep02_March.pdf",width=5, height=5)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 1, 0, 0))
boxplot(I(value-0.5)~censor + method, data = dfA, xaxt = "n", main = "MLE/CensMLE")
abline(h = 0, col = "red")
abline(v = 4.5, col = "grey")
axis(side=1, at= 1:8, labels= c("36", "60", "90", "120", "36", "60", "90", "120"))
mtext(text="Censoring Time", side=1, line=4, las=1)
mtext(text="Observed Bias", side=2, line=4, las=0)
dev.off()








