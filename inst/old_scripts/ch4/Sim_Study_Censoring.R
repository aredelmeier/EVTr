
# This R code calculates the mean injury length for non-censored injuries vs censored injuries

setwd("/Users/Redelmeier/Dropbox/Research/Simulations/Feb 9 New Code")

library(simsalapar)
library(parallel)
library(foreach)
library(doParallel)
library(QRM)
library(dplyr)

source("/Users/Redelmeier/Dropbox/Research/Simulations/Feb 9 New Code/Source_keep.R")

################################################################  
####################### All Observations #######################
################################################################ 


################################################################ 
######################### Theta = 1 ############################
################################################################

varList <- 
  varlist(
    n.sim = list(type = "N", expr = quote(N[sim]), value = 1000),
    n = list(type = "frozen", value = 1000),
    num_inj = list(type = "frozen", value = 30),
    xi = list(type = "frozen", value = 1),
    censor = list(type = "grid", expr = quote(censor), value = c(3, 6, 12, 36)),
    beta = list(type = "frozen",  value = 1),
    rate_exp = list(type = "frozen", value = 1),
    method = list(type = "inner", expr = quote(method), value = c("Norm_Value", "Cens_Value")) 
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
clusterEvalQ(cl,source("/Users/Redelmeier/Dropbox/Research/Simulations/Feb 9 New Code/Source_keep.R"))

#res <- doForeach(varList, cluster=cl, seed="seq", doOne=doOne,  sfile="Cens_keep2.rds", monitor=TRUE)
res <- doForeach(varList, cluster=cl, seed="seq", doOne=doOne,  sfile="Cens_keep_361236.rds", monitor=TRUE)
#res <- doForeach(varList, cluster=cl, seed="seq", doOne=doOne,  sfile="Cens_keep_xi1.rds", monitor=TRUE)
val <- getArray(res)
df <- array2df(val)

df_Normal36 <- df %>% dplyr::filter(method == "Norm_Value", censor == "36")
df_Censor36 <- df %>% dplyr::filter(method == "Cens_Value", censor == "36")


setwd("/Users/Redelmeier/Dropbox/Research/Presentations/To_Show_Plots")

pdf("Sim3_hist_all.pdf",width=5, height=5)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 5, 0, 0))
hist(df_Normal36$value, main = "theta = 1, non-censored obs, c = 36", xlab = "")
mtext(text="Average Injury Length", side=1, line=4, las=1)
dev.off()

pdf("Sim3_hist_all2.pdf",width=5, height=5)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 5, 0, 0))
hist(df_Censor36$value, breaks=seq(0,160,by=2), main = "theta = 1, censored obs, c = 36", xlab = "")
mtext(text="Average Injury Length", side=1, line=4, las=1)
dev.off()

#setwd("/Users/Redelmeier/Dropbox/Research/Presentations/To_Show_Plots")
setwd("/Users/Redelmeier/Dropbox/Research/Thesis")

pdf("Sim3_1.pdf",width=6, height=5)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 5, 0, 0))
boxplot(I(value)~method + censor, data = df, xaxt = "n", ylim = c(0,25))
abline(v=2.5,col="grey")
abline(v=4.5, col = "grey")
abline(v=6.5, col = "grey")
axis(side=1, at= 1:8, label = c("NC", "C", "NC", "C", "NC", "C", "NC", "C"))
#axis(side=1, at= c(1.5, 3.5, 5.5, 7.5), label = c("3", "6", "12", "36"))
mtext(text= paste("  tau = 3 ","    ","   tau = 6 ","      "," tau = 12   ","  ", "tau = 36"), 
      side=1, line=2.5, las=1)
mtext(text="End-of-Followup Time", side=1, line=4, las=1)
mtext(text="Average Injury Length", side=2, line=4, las=0)
dev.off()

pdf("Sim3_1_2.pdf",width=5, height=5)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 5, 0, 0))
boxplot(I(value)~method + censor, data = df, xaxt = "n", main = "theta = 1", ylim = c(0, 50))
abline(v=2.5,col="grey")
abline(v=4.5, col = "grey")
axis(side=1, at= 1:6, label = c("12", "12", "24", "24", "36", "36"))
mtext(text="Censoring Time", side=1, line=4, las=1)
mtext(text="Average Injury Length", side=2, line=4, las=0)
dev.off()


################################################################ 
######################### Theta = 1 April ######################
################################################################

varList <- 
  varlist(
    n.sim = list(type = "N", expr = quote(N[sim]), value = 1000),
    n = list(type = "frozen", value = 1000),
    num_inj = list(type = "frozen", value = 60),
    xi = list(type = "frozen", value = 1),
    censor = list(type = "grid", expr = quote(censor), value = c(36, 60, 90)),
    beta = list(type = "frozen",  value = 1),
    rate_exp = list(type = "frozen", value = 1),
    method = list(type = "inner", expr = quote(method), value = c("Norm_Value", "Cens_Value")) 
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
clusterEvalQ(cl,source("/Users/Redelmeier/Dropbox/Research/Simulations/Feb 9 New Code/Source_keep.R"))


res <- doForeach(varList, cluster=cl, seed="seq", doOne=doOne,  sfile="Cens_keep_60.rds", monitor=TRUE)
val <- getArray(res)
df <- array2df(val)

df_Normal36 <- df %>% dplyr::filter(method == "Norm_Value", censor == "36")
df_Censor36 <- df %>% dplyr::filter(method == "Cens_Value", censor == "36")


setwd("/Users/Redelmeier/Dropbox/Research/Presentations/To_Show_Plots")

pdf("Sim3_hist_all.pdf",width=5, height=5)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 5, 0, 0))
hist(df_Normal36$value, main = "theta = 1, non-censored obs, c = 36", xlab = "")
mtext(text="Average Injury Length", side=1, line=4, las=1)
dev.off()

pdf("Sim3_hist_all2.pdf",width=5, height=5)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 5, 0, 0))
hist(df_Censor36$value, breaks=seq(0,160,by=2), main = "theta = 1, censored obs, c = 36", xlab = "")
mtext(text="Average Injury Length", side=1, line=4, las=1)
dev.off()

#setwd("/Users/Redelmeier/Dropbox/Research/Presentations/To_Show_Plots")
setwd("/Users/Redelmeier/Dropbox/Research/Thesis")

pdf("Sim3_1.pdf",width=6, height=5)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 5, 0, 0))
boxplot(I(value)~method + censor, data = df, xaxt = "n")
abline(v=2.5,col="grey")
abline(v=4.5, col = "grey")
axis(side=1, at= 1:6, label = c("36", "36", "60", "60", "90", "90"))
mtext(text="End-of-Followup Time", side=1, line=4, las=1)
mtext(text="Average Injury Length", side=2, line=4, las=0)
dev.off()

pdf("Sim3_1_2.pdf",width=5, height=5)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 5, 0, 0))
boxplot(I(value)~method + censor, data = df, xaxt = "n", main = "theta = 1", ylim = c(0, 50))
abline(v=2.5,col="grey")
abline(v=4.5, col = "grey")
axis(side=1, at= 1:6, label = c("12", "12", "24", "24", "36", "36"))
mtext(text="Censoring Time", side=1, line=4, las=1)
mtext(text="Average Injury Length", side=2, line=4, las=0)
dev.off()


################################################################ 
######################### Theta = 0.2 ##########################
################################################################

setwd("/Users/Redelmeier/Dropbox/Research/Simulations/Feb 9 New Code")

varList2 <- 
  varlist(
    n.sim = list(type = "N", expr = quote(N[sim]), value = 1000),
    n = list(type = "frozen", value = 1000),
    num_inj = list(type = "frozen", value = 30),
    xi = list(type = "frozen", value = 0.5),
    censor = list(type = "grid", expr = quote(censor), value = c(36, 90, 120)),
    beta = list(type = "frozen",  value = 1),
    rate_exp = list(type = "frozen", value = 0.2),
    method = list(type = "inner", expr = quote(method), value = c("Norm_Value", "Cens_Value")) 
  )

# Parallel computing
no_cores <- detectCores()-1
cl <- makeCluster(no_cores)
clusterExport(cl,"df")
clusterEvalQ(cl,library(QRM))
clusterEvalQ(cl,source("/Users/Redelmeier/Dropbox/Research/Simulations/Feb 9 New Code/Source_keep.R"))

res2 <- doForeach(varList2, cluster=cl, seed="seq", doOne=doOne,  sfile="Cens_keep.rds",monitor=TRUE)
val2 <- getArray(res2)
df2 <- array2df(val2)

df2_Normal120 <- df2 %>% dplyr::filter(method == "Norm_Value", censor == "120")
df2_Censor120 <- df2 %>% dplyr::filter(method == "Cens_Value", censor == "120")


setwd("/Users/Redelmeier/Dropbox/Research/Presentations/To_Show_Plots")

pdf("Sim3_hist_all_2.pdf",width=5, height=5)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 5, 0, 0))
hist(df2_Normal120$value, main = "theta = 0.2, non-censored obs, c = 120", xlab = "")
mtext(text="Average Injury Length", side=1, line=4, las=1)
dev.off()

pdf("Sim3_hist_all_2_2.pdf",width=5, height=5)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 5, 0, 0))
hist(df2_Censor120$value, breaks=seq(0,350,by=2), xlim = c(0, 250), main = "theta = 0.2, censored obs, c = 120", xlab = "")
mtext(text="Average Injury Length", side=1, line=4, las=1)
dev.off()

setwd("/Users/Redelmeier/Dropbox/Research/Presentations/To_Show_Plots")

pdf("Sim3_2.pdf",width=5, height=5)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 5, 0, 0))
boxplot(I(value)~method + censor, data = df2,  xaxt = "n", main = "theta = 0.2")
abline(v=2.5,col="grey")
abline(v=4.5, col = "grey")
axis(side=1, at= 1:6, label = c("36", "36", "90", "90", "120", "120"))
mtext(text="Censoring Time", side=1, line=4, las=1)
mtext(text="Average Injury Length", side=2, line=4, las=0)
dev.off()

pdf("Sim3_2_2.pdf",width=5, height=5)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 5, 0, 0))
boxplot(I(value)~method + censor, data = df2,  xaxt = "n", main = "theta = 0.2", ylim = c(0, 50))
abline(v=2.5,col="grey")
abline(v=4.5, col = "grey")
axis(side=1, at= 1:6, label =  c("36", "36", "90", "90", "120", "120"))
mtext(text="Censoring Time", side=1, line=4, las=1)
mtext(text="Average Injury Length", side=2, line=4, las=0)
dev.off()

################################################################  
####################### Just Maxima ############################
################################################################ 


################################################################ 
######################### Theta = 1 ############################
################################################################

setwd("/Users/Redelmeier/Dropbox/Research/Simulations/Feb 9 New Code")

varList3 <- 
  varlist(
    n.sim = list(type = "N", expr = quote(N[sim]), value = 1000),
    n = list(type = "frozen", value = 1000),
    num_inj = list(type = "frozen", value = 30),
    xi = list(type = "frozen", value = 0.5),
    censor = list(type = "grid", expr = quote(censor), value = c(24, 36)),
    beta = list(type = "frozen",  value = 1),
    rate_exp = list(type = "frozen", value = 1),
    ne = list(type = "frozen", value = 50),
    method = list(type = "inner", expr = quote(method), value = c("Norm_Value_max", "Cens_Value")) 
  )

doOne2 <- function(n, xi, censor, beta, rate_exp, num_inj, method, ne){
  stopifnot(require(QRM), require(dplyr), sapply(list(n, xi, censor, rate_exp, num_inj, beta,ne), is.numeric))
  met <- rep(NA, length(method))
  names(met) <- method
  
  full <- gen_data(censor, xi, n, num_inj, beta, rate_exp)
  
  m <- maximum(full)
  exc <- excess(m, ne)
  
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
clusterEvalQ(cl,source("/Users/Redelmeier/Dropbox/Research/Simulations/Feb 9 New Code/Source_keep.R"))

res3 <- doForeach(varList3, cluster=cl, seed="seq", doOne=doOne2,  sfile="Cens_keep3.rds",monitor=TRUE)
val3 <- getArray(res3)
df3 <- array2df(val3)

df3_Normal36 <- df3 %>% dplyr::filter(method == "Norm_Value_max", censor == "36")
df3_Censor36 <- df3 %>% dplyr::filter(method == "Cens_Value", censor == "36")

setwd("/Users/Redelmeier/Dropbox/Research/Presentations/To_Show_Plots")

pdf("Sim3_hist1.pdf",width=5, height=5)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 5, 0, 0))
hist(df3_Normal36$value, main = "theta = 1, non-censored maxima, 
     c = 36",xlab = "")
mtext(text="Average Injury Length", side=1, line=4, las=1)
dev.off()

pdf("Sim3_hist2.pdf",width=5, height=5)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 5, 0, 0))
hist(df3_Censor36$value, breaks=seq(0,4000,by=20), main = "theta = 1, censored maxima, 
     c = 36",xlab = "")
mtext(text="Average Injury Length", side=1, line=4, las=1)
dev.off()

setwd("/Users/Redelmeier/Dropbox/Research/Presentations/To_Show_Plots")

pdf("Sim3_max_1.pdf",width=5, height=5)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 5, 0, 0))
boxplot(I(value-0.5)~method + censor, data = df3,xaxt = "n", main = "theta = 1 *")
abline(v=2.5,col="grey")
axis(side=1, at= 1:4, label = c("24", "24", "36", "36"))
mtext(text="Censoring Time", side=1, line=4, las=1)
mtext(text="Average Injury Length", side=2, line=4, las=0)
dev.off()

pdf("Sim3_max_1_2.pdf",width=5, height=5)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 5, 0, 0))
boxplot(I(value-0.5)~method + censor, data = df3, xaxt = "n", main = "theta = 1 *", ylim = c(0, 400))
abline(v=2.5,col="grey")
axis(side=1, at= 1:4, label = c("24", "24", "36", "36"))
mtext(text="Censoring Time", side=1, line=4, las=1)
mtext(text="Average Injury Length", side=2, line=4, las=0)
dev.off()

################################################################ 
######################### Theta = 0.2 ##########################
################################################################

setwd("/Users/Redelmeier/Dropbox/Research/Simulations/Feb 9 New Code")

varList4 <- 
  varlist(
    n.sim = list(type = "N", expr = quote(N[sim]), value = 1000),
    n = list(type = "frozen", value = 1000),
    num_inj = list(type = "frozen", value = 30),
    xi = list(type = "frozen", value = 0.5),
    censor = list(type = "grid", expr = quote(censor), value = c(90, 120)),
    beta = list(type = "frozen",  value = 1),
    rate_exp = list(type = "frozen", value = 0.2),
    ne = list(type = "frozen", value = 50),
    method = list(type = "inner", expr = quote(method), value = c("Norm_Value_max", "Cens_Value")) 
  )


# Parallel computing
no_cores <- detectCores()-1
cl <- makeCluster(no_cores)
clusterExport(cl,"df")
clusterEvalQ(cl,library(QRM))
clusterEvalQ(cl,source("/Users/Redelmeier/Dropbox/Research/Simulations/Feb 9 New Code/Source_keep.R"))

res4 <- doForeach(varList4, cluster=cl, seed="seq", doOne=doOne2,  sfile="Cens_keep4.rds",monitor=TRUE)
val4 <- getArray(res4)
df4 <- array2df(val4)


df4_Normal120 <- df4 %>% dplyr::filter(method == "Norm_Value_max", censor == "120")
df4_Censor120 <- df4 %>% dplyr::filter(method == "Cens_Value", censor == "120")

setwd("/Users/Redelmeier/Dropbox/Research/Presentations/To_Show_Plots")

pdf("Sim3_hist_2.pdf",width=5, height=5)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 5, 0, 0))
hist(df4_Normal120$value, main = "theta = 0.2, non-censored maxima,
     c = 120",xlab = "")
mtext(text="Average Injury Length", side=1, line=4, las=1)
dev.off()

pdf("Sim3_hist_2_2.pdf",width=5, height=5)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 5, 0, 0))
hist(df4_Censor120$value, breaks=seq(0,6000,by=50), main = "theta = 0.2, censored maxima, 
     c = 120",xlab = "")
mtext(text="Average Injury Length", side=1, line=4, las=1)
dev.off()


setwd("/Users/Redelmeier/Dropbox/Research/Presentations/To_Show_Plots")

pdf("Sim3_max_2.pdf",width=5, height=5)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 5, 0, 0))
boxplot(I(value-0.5)~method + censor, data = df4, main = "theta = 0.2 *", xaxt = "n")
abline(v=2.5,col="grey")
axis(side=1, at= 1:4, label = c( "90", "90", "120", "120"))
mtext(text="Censoring Time", side=1, line=4, las=1)
mtext(text="Average Injury Length", side=2, line=4, las=0)
dev.off()

pdf("Sim3_max_2_2.pdf",width=5, height=5)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 5, 0, 0))
boxplot(I(value-0.5)~method + censor, data = df4, main = "theta = 0.2", xaxt = "n", ylim = c(0, 400))
abline(v=2.5,col="grey")
axis(side=1, at= 1:4, label = c("90", "90", "120", "120"))
mtext(text="Censoring Time", side=1, line=4, las=1)
mtext(text="Average Injury Length", side=2, line=4, las=0)
dev.off()


