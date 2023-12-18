#Library-------------------------------------------------------------------
pacman::p_load(netstat, RSelenium, rio,rvest,magrittr,dtplyr,forecast,DT,zoo,lubridate,hrbrthemes, data.table,arrow, seasonal,seasonalview,readxl, tools, xts,glmnet,dplyr,POET,stcov,fastICA,nlshrink,rmgarch,moments, PerformanceAnalytics, quadprog,NMOF, riskParityPortfolio, tidyquant,timetk,ggplot2,foreach,MASS,CovTools,readr,paran,factoextra,tidyverse,scales,tidyr,corrplot,reshape2,gridExtra,palmerpenguins,SciViews,dfms,boot,PeerPerformance,ggthemes)


#Data Importation--------------------------------------------------------------------
setwd("J:/C_Piedboeuf/Doc Admin/ING/Business Panel/Data")



# URL of the Investing.com page
url <- "https://fr.investing.com/indices/bel-20-historical-data"

# Read the HTML content of the webpage
webpage <- read_html(url)

# Extract the historical data table using CSS selector
historical_data <- webpage %>%
  html_node("#curr_table") %>%
  html_table()



#Data source: https://fr.investing.com/indices/bel-20-historical-data
BEL20 <- fread("BEL 20 - Données Historiques (1).csv")[11:167,]
ABI   <- fread("ABI - Données Historiques.csv")[1:97,]
ACKB  <- fread("ACKB - Données Historiques.csv")[1:97,]
AGES  <- fread("AGES - Données Historiques.csv")[1:97,]
AOO   <- fread("AOO - Données Historiques.csv")[1:97,]
APAM  <- fread("APAM - Données Historiques.csv")[1:97,]
ARGX  <- fread("ARGX - Données Historiques.csv")[1:97,]
BAR   <- fread("BAR - Données Historiques.csv")[1:97,]
COFB  <- fread("COFB - Données Historiques.csv")[1:97,]
ELI   <- fread("ELI - Données Historiques.csv")[1:97,]
GBLB  <- fread("GBLB - Données Historiques.csv")[1:97,]
GLPG  <- fread("GLPG - Données Historiques.csv")[1:97,]
IETB  <- fread("IETB - Données Historiques.csv")[1:97,]
KBC   <- fread("KBC - Données Historiques.csv")[1:97,]
MLXS  <- fread("MLXS - Données Historiques.csv")[1:97,]
PROX  <- fread("PROX - Données Historiques.csv")[1:97,]
SOF   <- fread("SOF - Données Historiques.csv")[1:97,]
SOLB  <- fread("SOLB - Données Historiques.csv")[1:97,]
UCB   <- fread("UCB - Données Historiques.csv")[1:97,]
UMI   <- fread("UMI - Données Historiques.csv")[1:97,]
WDPP  <- fread("WDPP - Données Historiques.csv")[1:97,]

nrow(ABI)
nrow(ACKB)
nrow(AGES)
nrow(AOO)
nrow(APAM) #144
nrow(ARGX) #102
nrow(BAR)
nrow(ELI)
nrow(GBLB)
nrow(GLPG)
nrow(IETB)
nrow(KBC)
nrow(MLXS)
nrow(PROX)
nrow(SOF)
nrow(SOLB)
nrow(UCB)
nrow(UMI)
nrow(WDPP)

P10 = read_delim("10_Industry_Portfolios.csv",delim = ";", escape_double = FALSE, trim_ws = TRUE)
P25 = read_delim("25_Portfolios_5x5.csv",delim = ";", escape_double = FALSE, trim_ws = TRUE)
P48 = read_delim("48_Industry_Portfolios.csv",delim = ";", escape_double = FALSE, trim_ws = TRUE)
P100= read_delim("100_Portfolios_10x10.csv",delim = ";", escape_double = FALSE, trim_ws = TRUE)

P10 = as.matrix(sapply(P10[523:1159,2:11], as.numeric))/100
P48 = as.matrix(sapply(P48[523:1159,2:49], as.numeric))/100
P25 = as.matrix(sapply(P25[523:1159,2:26], as.numeric))/100
P100= as.matrix(sapply(P100[523:1159,2:101], as.numeric))/100

#########################################################################################
######################### Minimum-variance portfolio ####################################
#########################################################################################

#Weights
W_MV_10 = array(dim=c(10,86)) 
W_MV_48 = array(dim=c(48,86))
W_MV_25 = array(dim=c(25,86))
W_MV_100= array(dim=c(100,86))

a=0
b=0
for(i in 1:86)
{
  a[i]=i*6-5
  b[i]=114+i*6
  W_MV_10[,i] = minvar(cov(P10[a[i]:b[i],]))
  W_MV_48[,i] = minvar(cov(P48[a[i]:b[i],]))
  W_MV_25[,i] = minvar(cov(P25[a[i]:b[i],]))
  W_MV_100[,i]= minvar(cov(P100[a[i]:b[i],]))
  
}

#Returns
R_MV_10 = array(dim=c(6,86))
R_MV_48 = array(dim=c(6,86))
R_MV_25 = array(dim=c(6,86))
R_MV_100= array(dim=c(6,86))

for (i in 1:86){
  a[i]=115+i*6
  b[i]=120+i*6
  
  R_MV_10[,i]<-(P10[a[i]:b[i],]%*%W_MV_10[,i])
  R_MV_48[,i]<-(P48[a[i]:b[i],]%*%W_MV_48[,i])
  R_MV_25[,i]<-(P25[a[i]:b[i],]%*%W_MV_25[,i])
  R_MV_100[,i]<-(P100[a[i]:b[i],]%*%W_MV_100[,i])
  
}

AR_MV_10 <- mean(R_MV_10)*12
AR_MV_48 <- mean(R_MV_48)*12
AR_MV_25 <- mean(R_MV_25)*12
AR_MV_100<- mean(R_MV_100)*12

SD_MV_10<-sd(R_MV_10)*sqrt(12)
SD_MV_48<-sd(R_MV_48)*sqrt(12)
SD_MV_25<-sd(R_MV_25)*sqrt(12)
SD_MV_100<-sd(R_MV_100)*sqrt(12)

SR_MV_10 = AR_MV_10/SD_MV_10  #1.00383
SR_MV_48 = AR_MV_48/SD_MV_48  #0.9497646
SR_MV_25 = AR_MV_25/SD_MV_25  #0.8888929
SR_MV_100= AR_MV_100/SD_MV_100 #0.7626527




#########################################################################################
######################### Mean-variance portfolio #######################################
#########################################################################################

#Weights
W_MSR_10 = array(dim=c(10,86)) 
W_MSR_48 = array(dim=c(48,86))
W_MSR_25 = array(dim=c(25,86))
W_MSR_100= array(dim=c(100,86))

a=0
b=0
for(i in 1:86)
{
  a[i]=i*6-5
  b[i]=114+i*6
  
  W_MSR_10[,i] = (solve(cov(P10[a[i]:b[i],]))%*%colMeans(P10[a[i]:b[i],]))/as.numeric(t(array(1,dim=c(10,1)))%*%solve(cov(P10[a[i]:b[i],]))%*%colMeans(P10[a[i]:b[i],]))
  W_MSR_48[,i] = (solve(cov(P48[a[i]:b[i],]))%*%colMeans(P48[a[i]:b[i],]))/as.numeric(t(array(1,dim=c(48,1)))%*%solve(cov(P48[a[i]:b[i],]))%*%colMeans(P48[a[i]:b[i],]))
  W_MSR_25[,i] = (solve(cov(P25[a[i]:b[i],]))%*%colMeans(P25[a[i]:b[i],]))/as.numeric(t(array(1,dim=c(25,1)))%*%solve(cov(P25[a[i]:b[i],]))%*%colMeans(P25[a[i]:b[i],]))
  W_MSR_100[,i]= (solve(cov(P100[a[i]:b[i],]))%*%colMeans(P100[a[i]:b[i],]))/as.numeric(t(array(1,dim=c(100,1)))%*%solve(cov(P100[a[i]:b[i],]))%*%colMeans(P100[a[i]:b[i],]))
  
}

#Returns
R_MSR_10 = array(dim=c(6,86))
R_MSR_48 = array(dim=c(6,86))
R_MSR_25 = array(dim=c(6,86))
R_MSR_100= array(dim=c(6,86))


for (i in 1:86){
  a[i]=115+i*6
  b[i]=120+i*6
  
  R_MSR_10[,i] = (P10[a[i]:b[i],]%*%W_MSR_10[,i])
  R_MSR_48[,i] = (P48[a[i]:b[i],]%*%W_MSR_48[,i])
  R_MSR_25[,i] = (P25[a[i]:b[i],]%*%W_MSR_25[,i])
  R_MSR_100[,i]= (P100[a[i]:b[i],]%*%W_MSR_100[,i])
  
 
}
AR_MSR_10 = mean(R_MSR_10)*12
AR_MSR_48 = mean(R_MSR_48)*12
AR_MSR_25 = mean(R_MSR_25)*12
AR_MSR_100= mean(R_MSR_100)*12

SD_MSR_10 = sd(R_MSR_10)*sqrt(12)
SD_MSR_48 = sd(R_MSR_48)*sqrt(12)
SD_MSR_25 = sd(R_MSR_25)*sqrt(12)
SD_MSR_100= sd(R_MSR_100)*sqrt(12)

SR_MSR_10 = AR_MSR_10/SD_MSR_10   #0.6704312
SR_MSR_48 = AR_MSR_48/SD_MSR_48   #0.3716232
SR_MSR_25 = AR_MSR_25/SD_MSR_25   #1.182515
SR_MSR_100= AR_MSR_100/SD_MSR_100 #0.0983142



#########################################################################################
######################### Equally-weighted portfolio ####################################
#########################################################################################

##Out-of-Sample
#Weights
W_EW_10 = matrix(1/10,  dim(P10)[2] )
W_EW_48 = matrix(1/48,  dim(P48)[2] )
W_EW_25 = matrix(1/25,  dim(P25)[2] )
W_EW_100= matrix(1/100, dim(P100)[2])

#Returns
R_EW_10 = array(dim=c(6,86))
R_EW_48 = array(dim=c(6,86))
R_EW_25 = array(dim=c(6,86))
R_EW_100= array(dim=c(6,86))

a=0
b=0
for (i in 1:86){
  a[i]=115+i*6
  b[i]=120+i*6
  
  R_EW_10[,i] =(P10[a[i]:b[i],]%*%W_EW_10)
  R_EW_48[,i] =(P48[a[i]:b[i],]%*%W_EW_48)
  R_EW_25[,i] =(P25[a[i]:b[i],]%*%W_EW_25)
  R_EW_100[,i]=(P100[a[i]:b[i],]%*%W_EW_100)
  
}

AR_EW_10 = mean(R_EW_10)*12
AR_EW_48 = mean(R_EW_48)*12
AR_EW_25 = mean(R_EW_25)*12
AR_EW_100= mean(R_EW_100)*12

SD_EW_10 = sd(R_EW_10)*sqrt(12)
SD_EW_48 = sd(R_EW_48)*sqrt(12)
SD_EW_25 = sd(R_EW_25)*sqrt(12)
SD_EW_100= sd(R_EW_100)*sqrt(12)

SR_EW_10 = AR_EW_10/SD_EW_10   #0.8659715
SR_EW_48 = AR_EW_48/SD_EW_48   #0.7782867
SR_EW_25 = AR_EW_25/SD_EW_25   #0.7598897
SR_EW_100= AR_EW_100/SD_EW_100 #0.6089611


#########################################################################################
######################### Asset-Variance-Parity portfolio ###############################
#########################################################################################

#PACKAGE METHOD
W_AVP_10 = riskParityPortfolio(cov(P10))
W_AVP_10 = W_AVP_10$w
W_AVP_25 = riskParityPortfolio(cov(P25))
W_AVP_25 = W_AVP_25$w
W_AVP_48 = riskParityPortfolio(cov(P48))
W_AVP_48 = W_AVP_48$w
W_AVP_100= riskParityPortfolio(cov(P100))
W_AVP_100= W_AVP_100$w


#OPTIMISATION METHOD
w_function = function(cov_matrix)
{
# Define the objective function
obj_func <- function(w, cov_matrix) {
  return(w %*% cov_matrix %*% w - sum(log(w)))
}

# Define the initial value for w
n_assets <- ncol(cov_matrix)
w_init <- rep(1/n_assets, n_assets)

# Set the optimization control parameters
ctrl <- list(fnscale = 1,  # minimise the objective function
             factr = 1e-8, # absolute tolerance for convergence
             maxit = 1000) # maximum number of iterations

# Optimize the objective function subject to the constraint
result <- optim(w_init, obj_func, cov_matrix = cov_matrix, method = "L-BFGS-B", lower = rep(0, n_assets), control = ctrl)

# Extract the optimal solution
w_optimal <- result$par/sum(result$par)
return(w_optimal)



}
#look at the risk proportion of every asset
W_optimal = w_function(P10)
A =array(dim = c(1,10))
for (i in 1:10)
{
  A[i] =(w_optimal[i] %*% (cov_matrix %*% w_optimal)[i])/(w_optimal %*% cov_matrix %*% w_optimal)
  
  
}

#Weights
W_AVP_10 = array(dim=c(10,86)) 
W_AVP_48 = array(dim=c(48,86))
W_AVP_25 = array(dim=c(25,86))
W_AVP_100= array(dim=c(100,86))

a=0
b=0
for(i in 1:86)
{
  a[i]=i*6-5
  b[i]=114+i*6
  
  W_AVP_10[,i] = t(w_function(cov(P10[a[i]:b[i],])))
  W_AVP_48[,i] = t(w_function(cov(P48[a[i]:b[i],])))
  W_AVP_25[,i] = t(w_function(cov(P25[a[i]:b[i],])))
  W_AVP_100[,i]= t(w_function(cov(P100[a[i]:b[i],])))
  
}

#Returns
R_AVP_10 = array(dim=c(6,86))
R_AVP_48 = array(dim=c(6,86))
R_AVP_25 = array(dim=c(6,86))
R_AVP_100= array(dim=c(6,86))

for (i in 1:86){
  a[i]=115+i*6
  b[i]=120+i*6
  
  R_AVP_10[,i]<-(P10[a[i]:b[i],]%*%W_AVP_10[,i])
  R_AVP_48[,i]<-(P48[a[i]:b[i],]%*%W_AVP_48[,i])
  R_AVP_25[,i]<-(P25[a[i]:b[i],]%*%W_AVP_25[,i])
  R_AVP_100[,i]<-(P100[a[i]:b[i],]%*%W_AVP_100[,i])
  
  
}

AR_AVP_10 = mean(R_AVP_10)*12
AR_AVP_48 = mean(R_AVP_48)*12
AR_AVP_25 = mean(R_AVP_25)*12
AR_AVP_100= mean(R_AVP_100)*12

SD_AVP_10 = sd(R_AVP_10)*sqrt(12)
SD_AVP_48 = sd(R_AVP_48)*sqrt(12)
SD_AVP_25 = sd(R_AVP_25)*sqrt(12)
SD_AVP_100= sd(R_AVP_100)*sqrt(12)

SR_AVP_10 = AR_AVP_10/SD_AVP_10   #0.9148238
SR_AVP_48 = AR_AVP_48/SD_AVP_48   #0.8249968
SR_AVP_25 = AR_AVP_25/SD_AVP_25   #0.7821915
SR_AVP_100= AR_AVP_100/SD_AVP_100 #0.6668529

#########################################################################################
######################### Ledoit-Wolf Covariance Estimator ##############################
#########################################################################################

## LW mean-variance
#weights

W_LWMSR_10 = array(dim=c(10,86))
W_LWMSR_48 = array(dim=c(48,86))
W_LWMSR_25 = array(dim=c(25,86))
W_LWMSR_100= array(dim=c(100,86))

a=0
b=0
for(i in 1:86)
{
  a[i]=i*6-5
  b[i]=114+i*6
  
  
  "W_LWMSR_10[,i] = (solve(CovEst.2003LW((P10[a[i]:b[i],]))$S)%*%colMeans(P10[a[i]:b[i],]))/as.numeric(t(array(1,dim=c(10,1)))%*%solve(CovEst.2003LW((P10[a[i]:b[i],]))$S)%*%colMeans(P10[a[i]:b[i],]))
  W_LWMSR_48[,i] = (solve(CovEst.2003LW((P48[a[i]:b[i],]))$S)%*%colMeans(P48[a[i]:b[i],]))/as.numeric(t(array(1,dim=c(48,1)))%*%solve(CovEst.2003LW((P48[a[i]:b[i],]))$S)%*%colMeans(P48[a[i]:b[i],]))
  W_LWMSR_25[,i] = (solve(CovEst.2003LW((P25[a[i]:b[i],]))$S)%*%colMeans(P25[a[i]:b[i],]))/as.numeric(t(array(1,dim=c(25,1)))%*%solve(CovEst.2003LW((P25[a[i]:b[i],]))$S)%*%colMeans(P25[a[i]:b[i],]))
  W_LWMSR_100[,i]= (solve(CovEst.2003LW((P100[a[i]:b[i],]))$S)%*%colMeans(P100[a[i]:b[i],]))/as.numeric(t(array(1,dim=c(100,1)))%*%solve(CovEst.2003LW((P100[a[i]:b[i],]))$S)%*%colMeans(P100[a[i]:b[i],]))"
  
  W_LWMSR_10[,i] = (solve(linshrink_cov((P10[a[i]:b[i],])))%*%colMeans(P10[a[i]:b[i],]))/as.numeric(t(array(1,dim=c(10,1)))%*%solve(linshrink_cov((P10[a[i]:b[i],])))%*%colMeans(P10[a[i]:b[i],]))
  W_LWMSR_48[,i] = (solve(linshrink_cov((P48[a[i]:b[i],])))%*%colMeans(P48[a[i]:b[i],]))/as.numeric(t(array(1,dim=c(48,1)))%*%solve(linshrink_cov((P48[a[i]:b[i],])))%*%colMeans(P48[a[i]:b[i],]))
  W_LWMSR_25[,i] = (solve(linshrink_cov((P25[a[i]:b[i],])))%*%colMeans(P25[a[i]:b[i],]))/as.numeric(t(array(1,dim=c(25,1)))%*%solve(linshrink_cov((P25[a[i]:b[i],])))%*%colMeans(P25[a[i]:b[i],]))
  W_LWMSR_100[,i]= (solve(linshrink_cov((P100[a[i]:b[i],])))%*%colMeans(P100[a[i]:b[i],]))/as.numeric(t(array(1,dim=c(100,1)))%*%solve(linshrink_cov((P100[a[i]:b[i],])))%*%colMeans(P100[a[i]:b[i],]))
  
  
}

#Returns
R_LWMSR_10 = array(dim=c(6,86))
R_LWMSR_48 = array(dim=c(6,86))
R_LWMSR_25 = array(dim=c(6,86))
R_LWMSR_100= array(dim=c(6,86))


for (i in 1:86){
  a[i]=115+i*6
  b[i]=120+i*6
  
  R_LWMSR_10[,i] = (P10[a[i]:b[i],]%*%W_LWMSR_10[,i])
  R_LWMSR_48[,i] = (P48[a[i]:b[i],]%*%W_LWMSR_48[,i])
  R_LWMSR_25[,i] = (P25[a[i]:b[i],]%*%W_LWMSR_25[,i])
  R_LWMSR_100[,i]= (P100[a[i]:b[i],]%*%W_LWMSR_100[,i])
  
  
}
AR_LWMSR_10 = mean(R_LWMSR_10)*12
AR_LWMSR_48 = mean(R_LWMSR_48)*12
AR_LWMSR_25 = mean(R_LWMSR_25)*12
AR_LWMSR_100= mean(R_LWMSR_100)*12

SD_LWMSR_10 = sd(R_LWMSR_10)*sqrt(12)
SD_LWMSR_48 = sd(R_LWMSR_48)*sqrt(12)
SD_LWMSR_25 = sd(R_LWMSR_25)*sqrt(12)
SD_LWMSR_100= sd(R_LWMSR_100)*sqrt(12)

SR_LWMSR_10 = AR_LWMSR_10/SD_LWMSR_10   #0.7163915
SR_LWMSR_48 = AR_LWMSR_48/SD_LWMSR_48   #0.446826
SR_LWMSR_25 = AR_LWMSR_25/SD_LWMSR_25   #1.273495
SR_LWMSR_100= AR_LWMSR_100/SD_LWMSR_100 #1.720445

boxplot(t(W_LWMSR_10), main = "Ledoit-Wolf Mean-Variance - P10")
boxplot(t(W_LWMSR_48), main = "Ledoit-Wolf Mean-Variance - P48")
boxplot(t(W_LWMSR_25), main = "Ledoit-Wolf Mean-Variance - P25")
boxplot(t(W_LWMSR_100),main= "Ledoit-Wolf Mean-Variance - P100")

## LW minimum-variance
#weights

W_LWMV_10 = array(dim=c(10,86))
W_LWMV_48 = array(dim=c(48,86))
W_LWMV_25 = array(dim=c(25,86))
W_LWMV_100= array(dim=c(100,86))

a=0
b=0
for(i in 1:86)
{
  a[i]=i*6-5
  b[i]=114+i*6
  
  W_LWMV_10[,i] = minvar(linshrink_cov(((P10[a[i]:b[i],]))))
  W_LWMV_48[,i] = minvar(linshrink_cov(((P48[a[i]:b[i],]))))
  W_LWMV_25[,i] = minvar(linshrink_cov(((P25[a[i]:b[i],]))))
  W_LWMV_100[,i]= minvar(linshrink_cov(((P100[a[i]:b[i],]))))

  "W_LWMV_10[,i] = minvar(CovEst.2003LW(((P10[a[i]:b[i],])))$S)
  W_LWMV_48[,i] = minvar(CovEst.2003LW(((P48[a[i]:b[i],])))$S)
  W_LWMV_25[,i] = minvar(CovEst.2003LW(((P25[a[i]:b[i],])))$S)
  W_LWMV_100[,i]= minvar(CovEst.2003LW(((P100[a[i]:b[i],])))$S)"
}

#Returns
R_LWMV_10 = array(dim=c(6,86))
R_LWMV_48 = array(dim=c(6,86))
R_LWMV_25 = array(dim=c(6,86))
R_LWMV_100= array(dim=c(6,86))


for (i in 1:86){
  a[i]=115+i*6
  b[i]=120+i*6
  
  R_LWMV_10[,i] = (P10[a[i]:b[i],]%*%W_LWMV_10[,i])
  R_LWMV_48[,i] = (P48[a[i]:b[i],]%*%W_LWMV_48[,i])
  R_LWMV_25[,i] = (P25[a[i]:b[i],]%*%W_LWMV_25[,i])
  R_LWMV_100[,i]= (P100[a[i]:b[i],]%*%W_LWMV_100[,i])
  
  
}
AR_LWMV_10 = mean(R_LWMV_10)*12
AR_LWMV_48 = mean(R_LWMV_48)*12
AR_LWMV_25 = mean(R_LWMV_25)*12
AR_LWMV_100= mean(R_LWMV_100)*12

SD_LWMV_10 = sd(R_LWMV_10)*sqrt(12)
SD_LWMV_48 = sd(R_LWMV_48)*sqrt(12)
SD_LWMV_25 = sd(R_LWMV_25)*sqrt(12)
SD_LWMV_100= sd(R_LWMV_100)*sqrt(12)

SR_LWMV_10 = AR_LWMV_10/SD_LWMV_10   #1.011437
SR_LWMV_48 = AR_LWMV_48/SD_LWMV_48   #0.9585021
SR_LWMV_25 = AR_LWMV_25/SD_LWMV_25   #0.8874716
SR_LWMV_100= AR_LWMV_100/SD_LWMV_100 #0.773922

boxplot(t(W_LWMV_10), main = "Ledoit-Wolf Mean-Variance - P10",ylim = c(0, 0.5), col = "pink1")
boxplot(t(W_LWMV_48), main = "Ledoit-Wolf Mean-Variance - P48",ylim = c(0, 0.5))
boxplot(t(W_LWMV_25), main = "Ledoit-Wolf Mean-Variance - P25",ylim = c(0, 0.5))
boxplot(t(W_LWMV_100), main= "Ledoit-Wolf Mean-Variance - P100",ylim = c(0, 0.5))


#########################################################################################
######################### Optimal number of PCs (k) #####################################
#########################################################################################


##Bai & Ng -1-
ICr(P10)
ic10 = 9
ICr(P48)
ic48 = 4
ICr(P25)
ic25 = 20
ICr(P100)
ic100=6

##Bai & Ng -2-
baingcriterion <- function(X, rmax, rvar = NULL)
{
  X <- as.matrix(X)
  if(is.null(rvar)) rvar = rmax
  nobs = dim(X)[1]
  nvar = dim(X)[2]
  #mu <- matrix(1, nrow=nobs,ncol=1)%*%apply(X, 2, mean)
  #Winv <- diag(apply(X, 2, stats::sd)^-1)
  #x = (X-mu)%*%Winv
  x =scale(X)
  
  Gamma0 = stats::cov(x)
  H <- eigen(Gamma0)$vectors
  V <- vector()
  for(j in 1:rmax){
    I = x%*%(diag(nvar) - tcrossprod(H[,1:j]))
    V[j] = sum(I^2)/(nvar*nobs)
  }
  e = (nvar + nobs)/(nvar*nobs)
  
  penalty1 = (1:rmax)*e*log(1/e)
  penalty2 = (1:rmax)*e*log(min(nvar,nobs))
  penalty3 = (1:rmax)*(log(min(nvar,nobs))/min(nvar,nobs))
  
  PC = cbind(log(V) + penalty1, log(V) + penalty2, log(V) + penalty3)
  PC.min <- apply(PC, 2, which.min)
  return(list(PC = PC.min))
}

baingcriterion(PTest10Daily,10)


##Bai & Ng -3-
K_function = function(Data)
{
  Bound = 1000000000
  K_optimal = 0 
  Eig = eigen(cov(Data))
  V = Eig$vectors
  T = nrow(Data)
  N = ncol(Data)
  
  for( K in ncol(Data))
  {
  e = (t(Data) - V[,1:K]%*%t(V[,1:K])%*%(t(Data)))
  
  K_Op = ln(sum((e^2)/(N*T))) + K*((N+T)/(N*T))*ln((N*T)/(N+T))
  K_Op
  
  if ( K_Op < Bound)
  {
    Bound = K_Op
    K_optimal = K
  }
  K_optimal
  }
  
   return(K_optimal)
}



## Percentage of explained variance (>0.90)
K_EV_10 = array(dim=c(10,86))
K_EV_48 = array(dim=c(48,86))
K_EV_25 = array(dim=c(25,86))
K_EV_100= array(dim=c(100,86))

a=0
b=0

for(i in 1:86)
{
  a[i]=i*6-5
  b[i]=114+i*6
  
  eig_10 = eigen(cov(scale(P10[a[i]:b[i],])))
  eig_48 = eigen(cov(scale(P48[a[i]:b[i],])))
  eig_25 = eigen(cov(scale(P25[a[i]:b[i],])))
  eig_100= eigen(cov(scale(P100[a[i]:b[i],])))
  
  for (j in 1:10) {
    K_EV_10[j,i] = sum(eig_10$values[1:j])/sum(eig_10$values)
      }
  for (j in 1:48) {
    K_EV_48[j,i] = sum(eig_48$values[1:j])/sum(eig_48$values)
  }
  for (j in 1:25) {
    K_EV_25[j,i] = sum(eig_25$values[1:j])/sum(eig_25$values)
    
  }
  for (j in 1:100) {
    K_EV_100[j,i] = sum(eig_100$values[1:j])/sum(eig_100$values)
  }
}
k10= array(dim=c(86))
k48= array(dim=c(86))
k25= array(dim=c(86))
k100= array(dim=c(86))
for (j in 1:86){
for (i in 1:10)
{
  if (K_EV_10[i,j]>0.90)
    {
    k10[j] = i
    break
    }
}
  }
for (j in 1:86) {
for (i in 1:48)
{
  if (K_EV_48[i,j]>0.90)
  {
    k48[j] = i
    break
  }
}
}
for (j in 1:86){
for (i in 1:25)
{
  if (K_EV_25[i,j]>0.90)
  {
    k25[j] = i
    break
  }
}
}
for (j in 1:86) {
for (i in 1:100)
{
  if (K_EV_100[i,j]>0.90)
  {
    k100[j] = i
    break
  }
}
}

#Plot of explained variance percentage of individual components 
par(mfrow=c(2,2))
fviz_eig(prcomp(P10), addlabels = TRUE, linecolor = "Red", ylim = c(0, 85), xlab = "PC's", main = "P10")
fviz_eig(prcomp(P48), addlabels = TRUE, linecolor = "Green", ylim = c(0,85 ), xlab = "PC's", main = "P48")
fviz_eig(prcomp(P25), addlabels = TRUE, linecolor = "Blue", ylim = c(0, 85), xlab = "C's", main = "P25")
fviz_eig(prcomp(P100), addlabels = TRUE, linecolor = "Black", ylim = c(0, 85), xlab = "PC's", main = "P100")


## Scree Plot - Parallel Analysis
par(mfrow=c(2,2))
paran(P10, iterations = 5000, centile = 0, quietly = FALSE, 
      status = TRUE, all = TRUE, cfa = TRUE, graph = TRUE, color = TRUE, 
      col = c("Blue"), lty = c(1, 2, 3), lwd = 1, legend = FALSE, 
      file = "", width = 640, height = 640, grdevice = "png", seed = 0)

paran(P48, iterations = 5000, centile = 0, quietly = FALSE, 
      status = TRUE, all = TRUE, cfa = TRUE, graph = TRUE, color = TRUE, 
      col = c("Blue"), lty = c(1, 2, 3), lwd = 1, legend = FALSE, 
      file = "", width = 640, height = 640, grdevice = "png", seed = 0)
paran(P25, iterations = 5000, centile = 0, quietly = FALSE, 
      status = TRUE, all = TRUE, cfa = TRUE, graph = TRUE, color = TRUE, 
      col = c("Blue"), lty = c(1, 2, 3), lwd = 1, legend = FALSE, 
      file = "", width = 640, height = 640, grdevice = "png", seed = 0)
paran(P100, iterations = 5000, centile = 0, quietly = FALSE, 
      status = TRUE, all = TRUE, cfa = TRUE, graph = TRUE, color = TRUE, 
      col = c("Blue"), lty = c(1, 2, 3), lwd = 1, legend = FALSE, 
      file = "", width = 640, height = 640, grdevice = "png", seed = 0)



#########################################################################################
######################### PC-Variance-Parity portfolio ##################################
#########################################################################################
PCVP = function(dataset, k){
  
  #different parameters 
  sigma = cov(dataset) #matrice cov 
  V = eigen(sigma)$vectors[,1:k] #eigenvectors 
  lambda = diag(eigen(sigma)$values[1:k]) #eigenvalues 
  #View(lambda)
  LAMBDA = solve(sqrt(lambda)) #lambda^-1/2
  
  LAMBDA <- as.matrix(LAMBDA)
  V <- as.matrix(V)
  dataset <- as.matrix(dataset)
  
  
  #Weights computation
  one = as.matrix(rep(1, length(dataset[1,])))
  sign_1 = sign(LAMBDA %*% t(V) %*% one)
  
  w = ( V %*% LAMBDA  %*% sign_1 ) / as.numeric( t(one) %*% V %*% LAMBDA  %*% sign_1)
  w   
  
} 


# Weights
W_PCVP_10 = array(dim=c(10,86))
W_PCVP_48 = array(dim=c(48,86))
W_PCVP_25 = array(dim=c(25,86))
W_PCVP_100= array(dim=c(100,86))

a=0
b=0
for(i in 1:86)
{
  a[i]=i*6-5
  b[i]=114+i*6
  
  #eig_10 = eigen(cov((P10[a[i]:b[i],])))
  #eig_48 = eigen(cov((P48[a[i]:b[i],])))
  #eig_25 = eigen(cov((P25[a[i]:b[i],])))
  #eig_100= eigen(cov((P100[a[i]:b[i],])))
  
  
  #One10 = sign((solve(diag(eig_10$values))^(1/2))%*%t(eig_10$vectors)%*%rep(1,10))
  #One48 = sign(((solve(diag(eig_48$values)))^(1/2))%*%t(eig_48$vectors)%*%rep(1,48))
  #One25 = sign(((solve(diag(eig_25$values)))^(1/2))%*%t(eig_25$vectors)%*%rep(1,25))
  #One100= sign(((solve(diag(eig_100$values)))^(1/2))%*%t(eig_100$vectors)%*%rep(1,100))
 
  
  
  #W_PCVP_10[,i] = as.numeric(solve(rep(1,10)%*%eig_10$vectors %*% (solve(diag(eig_10$values))^(1/2)) %*% One10))*(eig_10$vectors %*% ((solve(diag(eig_10$values))^1/2)) %*% One10)
  #W_PCVP_48[,i] = as.numeric(solve(rep(1,48)%*%eig_48$vectors %*% (solve(diag(eig_48$values))^(1/2)) %*% One48))*(eig_48$vectors %*% ((solve(diag(eig_48$values))^1/2)) %*% One48)
  #W_PCVP_25[,i] = as.numeric(solve(rep(1,25)%*%eig_25$vectors %*% (solve(diag(eig_25$values))^(1/2)) %*% One25))*(eig_25$vectors %*% ((solve(diag(eig_25$values))^1/2)) %*% One25)
  #W_PCVP_100[,i]= as.numeric(solve(rep(1,100)%*%eig_100$vectors %*% (solve(diag(eig_100$values))^(1/2)) %*% One100))*(eig_100$vectors %*% ((solve(diag(eig_100$values))^1/2)) %*% One100)
  
  W_PCVP_10[,i]= PCVP(P10[a[i]:b[i],],k10[i])
  W_PCVP_48[,i]= PCVP(P48[a[i]:b[i],],k48[i])
  W_PCVP_25[,i]= PCVP(P25[a[i]:b[i],],k25[i])
  W_PCVP_100[,i]= PCVP(P100[a[i]:b[i],],k100[i])
  
  
}

#Returns
R_PCVP_10 = array(dim=c(6,86))
R_PCVP_48 = array(dim=c(6,86))
R_PCVP_25 = array(dim=c(6,86))
R_PCVP_100= array(dim=c(6,86))

for (i in 1:86){
  a[i]=115+i*6
  b[i]=120+i*6
  
  R_PCVP_10[,i]<-(P10[a[i]:b[i],]%*%W_PCVP_10[,i])
  R_PCVP_48[,i]<-(P48[a[i]:b[i],]%*%W_PCVP_48[,i])
  R_PCVP_25[,i]<-(P25[a[i]:b[i],]%*%W_PCVP_25[,i])
  R_PCVP_100[,i]<-(P100[a[i]:b[i],]%*%W_PCVP_100[,i])
  
  
}

AR_PCVP_10 = mean(R_PCVP_10)*12      #0.1105249
AR_PCVP_48 = mean(R_PCVP_48)*12      #0.1562515
AR_PCVP_25 = mean(R_PCVP_25)*12      #0.1330714
AR_PCVP_100= mean(R_PCVP_100)*12     #0.1147836

SD_PCVP_10 = sd(R_PCVP_10)*sqrt(12)  #0.1471289
SD_PCVP_48 = sd(R_PCVP_48)*sqrt(12)  #0.1688184
SD_PCVP_25 = sd(R_PCVP_25)*sqrt(12)  #0.2073243
SD_PCVP_100= sd(R_PCVP_100)*sqrt(12) #0.1831451

SR_PCVP_10 = AR_PCVP_10/SD_PCVP_10   #0.7512114
SR_PCVP_48 = AR_PCVP_48/SD_PCVP_48   #0.9255593
SR_PCVP_25 = AR_PCVP_25/SD_PCVP_25   #0.6418515
SR_PCVP_100= AR_PCVP_100/SD_PCVP_100 #0.6267357

par(mfrow=c(2,2))
boxplot(t(W_PCVP_10), xlab ="Assets P10", ylab = "Weights", col = "bisque1")
boxplot(t(W_PCVP_25)[,1:10], xlab ="Assets P25", ylab = "Weights", col = "bisque1")
boxplot(t(W_PCVP_48)[,1:10], xlab ="Assets P48", ylab = "Weights", col = "bisque1")
boxplot(t(W_PCVP_100)[,1:10],xlab ="Assets P100", ylab = "Weights", col = "bisque1")


#########################################################################################
######################### PCA-Subspace portfolios #######################################
#########################################################################################

##PC-mean-variance
#weights

PC_MSR_function = function(Data, k)
{
W_PCMSR = array(dim=c(ncol(Data),86))

a=0
b=0
for(i in 1:86)
{
  a[i]=i*6-5
  b[i]=114+i*6
  eig = eigen(cov((Data[a[i]:b[i],])))
  W_PCMSR[,i] = as.numeric(solve(rep(1,ncol(Data))%*%eig$vectors[,1:k[i]] %*% solve(diag(eig$values[1:k[i]])) %*% t(eig$vectors[,1:k[i]])%*%colMeans(Data[a[i]:b[i],])))*(eig$vectors[,1:k[i]] %*% solve(diag(eig$values[1:k[i]])) %*% t(eig$vectors[,1:k[i]]) %*% colMeans(Data[a[i]:b[i],]))

  }

#boxplot(t(W_PCMSR)[,1:10],xlab ="Assets", ylab = "Weights", col = "bisque1")

#Returns
R_PCMSR = array(dim=c(6,86))
for (i in 1:86){
  a[i]=115+i*6
  b[i]=120+i*6
  R_PCMSR[,i] = (Data[a[i]:b[i],]%*%W_PCMSR[,i])
}
AR_PCMSR = mean(R_PCMSR)*12
SD_PCMSR = sd(R_PCMSR)*sqrt(12)
SR_PCMSR = AR_PCMSR/SD_PCMSR

return(SR_PCMSR)

}
par(mfrow=c(2,2))
PC_MSR_function(P10, k10)   #0.6599168
PC_MSR_function(P48, k48)   #0.4427141
PC_MSR_function(P25, k25)   #0.6898977
PC_MSR_function(P100, k100) #1.097311

#Change in function of k
par(mfrow=c(2,2))
seq = rep(1,86)
SR = array(dim=c(10))
for(i in 2:10){
SR[i] = PC_MSR_function(P10, i*seq)   
}
plot(SR, type = 'l', col = "lightblue", xlab = "Number of components", ylab = "SR", main = "PC-MSR (P10)")

SR = array(dim=c(48))
for(i in 2:48){
  SR[i] = PC_MSR_function(P48, i*seq)   
}
plot(SR, type = 'l', col = "lightblue", xlab = "Number of components", ylab = "SR", main = "PC-MSR (P48)")

SR = array(dim=c(25))
for(i in 2:25){
  SR[i] = PC_MSR_function(P25, i*seq)   
}
plot(SR, type = 'l', col = "lightblue", xlab = "Number of components", ylab = "SR", main = "PC-MSR (P25)")

SR = array(dim=c(100))
for(i in 2:100){
  SR[i] = PC_MSR_function(P100, i*seq)   
}
plot(SR, type = 'l', col = "lightblue", xlab = "Number of components", ylab = "SR", main = "PC-MSR (P100)")


##PC-minimum-variance
#weights

PC_MV_function = function(Data, k)
{
W_PCMV = array(dim=c(ncol(Data),86))

a=0
b=0
for(i in 1:86)
{
  a[i]=i*6-5
  b[i]=114+i*6
  eig = eigen(cov((Data[a[i]:b[i],])))
  W_PCMV[,i] = as.numeric(solve(rep(1,ncol(Data))%*%eig$vectors[,1:k[i]] %*% solve(diag(eig$values[1:k[i]])) %*% t(eig$vectors[,1:k[i]])%*%rep(1,ncol(Data))))*(eig$vectors[,1:k[i]] %*% solve(diag(eig$values[1:k[i]])) %*% t(eig$vectors[,1:k[i]]) %*% rep(1,ncol(Data)))
}

#boxplot(t(W_PCMV)[,1:10],xlab ="Assets", ylab = "Weights", col = "bisque1")
#Returns
R_PCMV = array(dim=c(6,86))

for (i in 1:86){
  a[i]=115+i*6
  b[i]=120+i*6
  
  R_PCMV[,i] = (Data[a[i]:b[i],]%*%W_PCMV[,i])
}

AR_PCMV = mean(R_PCMV)*12
SD_PCMV = sd(R_PCMV)*sqrt(12)
SR_PCMV = AR_PCMV/SD_PCMV  
return(AR_PCMV)
}
PC_MV_function(P10, k10)   #0.9841925
PC_MV_function(P48, k48)   #1.088623
PC_MV_function(P25, k25)   #0.728117
PC_MV_function(P100, k100) #1.008152

#Change in function of k
par(mfrow=c(2,2))
SR = array(dim=c(10))
seq = rep(1,86)
for(i in 2:10){
  SR[i] = PC_MV_function(P10, i*seq)   
}
plot(SR, type = 'l', col = "lightblue", xlab = "Number of components", ylab = "SR", main = "PC-MV (P10)")

SR = array(dim=c(48))
for(i in 2:48){
  SR[i] = PC_MV_function(P48, i*seq)   
}
plot(SR, type = 'l', col = "lightblue", xlab = "Number of components", ylab = "SR", main = "PC-MV (P48)")

SR = array(dim=c(25))
for(i in 2:25){
  SR[i] = PC_MV_function(P25, i*seq)   
}
plot(SR, type = 'l', col = "lightblue", xlab = "Number of components", ylab = "SR", main = "PC-MV (P25)")

SR = array(dim=c(100))
for(i in 2:100){
  SR[i] = PC_MV_function(P100, i*seq)   
}
plot(SR, type = 'l', col = "lightblue", xlab = "Number of components", ylab = "SR", main = "PC-MV (P100)")



#########################################################################################
######################### Severini-PC portfolio #########################################
#########################################################################################
#weights

SR_S_function = function(P10,neig){
a=0
b=0
W_S_10 = array(dim=c(ncol(P10),86))
for(i in 1:86)
{
  a[i]=i*6-5
  b[i]=114+i*6
  
  eig_10 = eigen(cov((P10[a[i]:b[i],])))
  
  W_S_10[,i] = (eig_10$vectors[,neig])/sum(eig_10$vectors[,neig])
  
    
 
  
}
#boxplot(t(W_S_10)[,1:10])

#Returns
R_S_10 = array(dim=c(6,86))

for (i in 1:86){
  a[i]=115+i*6
  b[i]=120+i*6
  
  R_S_10[,i]<-(P10[a[i]:b[i],]%*%W_S_10[,i])
 
  
  
}

AR_S_10 = mean(R_S_10)*12
SD_S_10 = sd(R_S_10)*sqrt(12)
SR_S_10 = AR_S_10/SD_S_10   
return(SR_S_10)
}

#This is quasi the same as the 1/N
SR_S_function(P10,1)  #0.7970284
SR_S_function(P48,1)  #0.7303696
SR_S_function(P25,1)  #0.7258577
SR_S_function(P100,1) #0.01331437

SR_S_function(P10,2)  #-0.2340491
SR_S_function(P48,2)  #0.1180155
SR_S_function(P25,2)  #0.1650339
SR_S_function(P100,2) #-0.1728867


S_efunction = function(Data){
SR_S_e = array(dim=c(ncol(Data)))
a=0
e10 = 0
for( i in 1:ncol(Data)){
  
  SR_S_e[i] = SR_S_function(Data,i)
  
}
return(SR_S_e)
}

S_eOptifunction = function(Data){
  SR_S_e = array(dim=c(ncol(Data)))
  a=0
  e10 = 0
  for( i in 1:ncol(Data)){
    
    SR_S_e[i] = SR_S_function(Data,i)
    if(a<SR_S_e[i])
    { 
      a= SR_S_e[i]
      e10 = i
      
    }
  }
  return(e10)
}
S_eOptifunction(P48)
S_eOptifunction(P10)
S_eOptifunction(P100)
S_eOptifunction(P25)

par(mfrow=c(2,2))
plot(S_efunction(P10), type = "l", xlab = "Number of eigenvectors", ylab = "SR", main = "Severini (P10)")
plot(S_efunction(P48), type = "l", xlab = "Number of eigenvectors", ylab = "SR", main = "Severini (P48)")
plot(S_efunction(P25), type = "l", xlab = "Number of eigenvectors", ylab = "SR", main = "Severini (P25)")
plot(S_efunction(P100), type = "l", xlab = "Number of eigenvectors", ylab = "SR", main = "Severini (P100)")

#########################################################################################
######################### Bounded-Noise portfolio #######################################
#########################################################################################

# Functions
covBN = function(data, gamma = 0.25, type = "CV", n_folds = 3, n_boot = 100){
  if(type == "CV"){
    idx_m = getIdxMFromCV(data, gamma, n_folds)
  }else{
    idx_m = getIdxMFromBoot(data, gamma, n_boot)
  }
  
  covNew(cov(data), idx_m$idx, idx_m$m)
}

covNew = function(cov_est, idx, m){
  n_stocks = ncol(cov_est)
  eigen_est = eigen(cov_est)
  eigen_est$values[idx:n_stocks] = eigen_est$values[idx:n_stocks] + m
  
  eigen_est$vectors %*% diag(eigen_est$values) %*% t(eigen_est$vectors)
}

# 3-fold cross validation (CV) is much quicker than bootstrap but less stable.
# Empirically speaking, the final portfolios (CV or bootstrap) has similar out-of-sample standard deviation.
getIdxMFromCV = function(data, gamma = 0.25, n_folds = 3){
  n_obs = nrow(data)
  n_stocks = ncol(data)
  
  group_size = ceiling(n_obs/n_folds)
  end_idx = min(ceiling(n_obs/n_folds * (n_folds - 1)) - 1, n_stocks)
  
  ratio_mat = matrix(NA, nrow = n_folds, ncol = end_idx)
  
  cov_est = cov(data)
  for(i in 1:n_folds){
    hold_out_idx = 1:group_size + (i-1) * group_size
    cov_in = cov(data[-hold_out_idx,])
    ratio_mat[i,] = getRealVarEstVarRatio(cov_est, cov_in, end_idx)
  }
  
  ratio_median = apply(ratio_mat, 2, median)
  
  noise_idx = which(ratio_median > 1 + gamma)
  
  if(length(noise_idx) == 0){
    idx = n_stocks
    m = 0
  }else{
    idx = min(noise_idx)
    m_vec = rep(NA, n_folds)
    
    for(i in 1:n_folds){
      hold_out_idx = 1:group_size + (i-1) * group_size
      cov_in = cov(data[-hold_out_idx,])
      
      eigen_in = eigen(cov_in)
      
      proj_matrix = t(eigen_in$vectors[,idx:n_stocks]) %*% (cov_est - cov_in) %*% eigen_in$vectors[,idx:n_stocks]
      
      m_vec[i] = eigen(proj_matrix)$values[1]
    }
    
    m = median(m_vec)
  }
  
  list(idx = idx, m = m)
}

getRealVarEstVarRatio = function(cov_true, cov_est, end_idx){
  eigen_est = eigen(cov_est)
  ratio = diag(t(eigen_est$vectors) %*% cov_true %*% eigen_est$vectors)/eigen_est$values
  ratio[1:end_idx]
}

# Weights
BNP_function = function(Data)
{
W_BNP = array(dim=c(ncol(Data),86))
a=0
b=0
for(i in 1:86)
{
  a[i]=i*6-5
  b[i]=114+i*6
  
  W_BNP[,i]  = rowSums(solve(covBN(Data[a[i]:b[i],], gamma = 0.25, type = "CV")))/sum(solve(covBN(Data[a[i]:b[i],], gamma = 0.25, type = "CV")))
}

#boxplot(t(W_BNP)[,1:10])
#Returns
R_BNP = array(dim=c(6,86))
i = 1

for (i in 1:86){
  a[i]=115+i*6
  b[i]=120+i*6
  
  R_BNP[,i]<-(Data[a[i]:b[i],]%*%W_BNP[,i])
  
}

AR = mean(R_BNP)*12
SD = sd(R_BNP)*sqrt(12)
SR = AR/SD
return(SR)
}

BNP_function(P10)  #0.9981173
BNP_function(P48)  #1.078119
BNP_function(P25)  #1.275393
BNP_function(P100) #1.318981


#########################################################################################
########################### Sharpe Ratio Test ###########################################
#########################################################################################



SR_Test_Function = function(R1, R2)
{
BiReturns = array(dim = c(86*6,2))
a = 0

for (i in 1:86)
{
  for ( j in 1:6)
  {
    a = a + 1
    BiReturns[a,1] = R1[j,i]
    BiReturns[a,2] = R2[j,i]
    
  }
  
}

return(block.size.calibrate(BiReturns))
}

#Change in the functions what it returns so that it returns the vector of return, this can now be used in the SR test
SR_Test_Function(R_PCVP_10,PC_MSR_function(P10,k10))
SR_Test_Function(R_PCVP_48,PC_MSR_function(P48,k48))
SR_Test_Function(R_PCVP_25,PC_MSR_function(P25,k25))
SR_Test_Function(R_PCVP_100,PC_MSR_function(P100,k100))

SR_Test_Function(R_PCVP_10, SR_S_function(P10,1))
SR_Test_Function(R_PCVP_48, SR_S_function(P48,1))
SR_Test_Function(R_PCVP_25, SR_S_function(P25,1))
SR_Test_Function(R_PCVP_100, SR_S_function(P100,1))

SR_Test_Function(R_PCVP_10, SR_S_function(P10,2))
SR_Test_Function(R_PCVP_48, SR_S_function(P48,2))
SR_Test_Function(R_PCVP_25, SR_S_function(P25,2))
SR_Test_Function(R_PCVP_100, SR_S_function(P100,2))

SR_Test_Function(R_PCVP_10, BNP_function(P10))
SR_Test_Function(R_PCVP_48, BNP_function(P48))
SR_Test_Function(R_PCVP_25, BNP_function(P25))
SR_Test_Function(R_PCVP_100, BNP_function(P100))

SR_Test_Function(R_PCVP_10,PC_MV_function(P10,k10))
SR_Test_Function(R_PCVP_48,PC_MV_function(P48,k48))
SR_Test_Function(R_PCVP_25,PC_MV_function(P25,k25))
SR_Test_Function(R_PCVP_100,PC_MV_function(P100,k100))

SR_Test_Function(PC_MV_function(P10,k10), PC_MSR_function(P10,k10))
SR_Test_Function(PC_MV_function(P48,k48), PC_MSR_function(P48,k48))
SR_Test_Function(PC_MV_function(P25,k25), PC_MSR_function(P25,k25))
SR_Test_Function(PC_MV_function(P100,k100), PC_MSR_function(P100,k100))

SR_Test_Function(PC_MV_function(P10,k10), SR_S_function(P10,1))
SR_Test_Function(PC_MV_function(P48,k48), SR_S_function(P48,1))
SR_Test_Function(PC_MV_function(P25,k25), SR_S_function(P25,1))
SR_Test_Function(PC_MV_function(P100,k100), SR_S_function(P100,1))

SR_Test_Function(PC_MV_function(P10,k10), SR_S_function(P10,2))
SR_Test_Function(PC_MV_function(P48,k48), SR_S_function(P48,2))
SR_Test_Function(PC_MV_function(P25,k25), SR_S_function(P25,2))
SR_Test_Function(PC_MV_function(P100,k100), SR_S_function(P100,2))

SR_Test_Function(PC_MV_function(P10,k10), BNP_function(P10))
SR_Test_Function(PC_MV_function(P48,k48), BNP_function(P48))
SR_Test_Function(PC_MV_function(P25,k25), BNP_function(P25))
SR_Test_Function(PC_MV_function(P100,k100), BNP_function(P100))

SR_Test_Function(PC_MSR_function(P10,k10), SR_S_function(P10,1))
SR_Test_Function(PC_MSR_function(P48,k48), SR_S_function(P48,1))
SR_Test_Function(PC_MSR_function(P25,k25), SR_S_function(P25,1))
SR_Test_Function(PC_MSR_function(P100,k100), SR_S_function(P100,1))

SR_Test_Function(PC_MSR_function(P10,k10), SR_S_function(P10,2))
SR_Test_Function(PC_MSR_function(P48,k48), SR_S_function(P48,2))
SR_Test_Function(PC_MSR_function(P25,k25), SR_S_function(P25,2))
SR_Test_Function(PC_MSR_function(P100,k100), SR_S_function(P100,2))

SR_Test_Function(PC_MSR_function(P10,k10), BNP_function(P10))
SR_Test_Function(PC_MSR_function(P48,k48), BNP_function(P48))
SR_Test_Function(PC_MSR_function(P25,k25), BNP_function(P25))
SR_Test_Function(PC_MSR_function(P100,k100), BNP_function(P100))

SR_Test_Function(SR_S_function(P10,1),BNP_function(P10))
SR_Test_Function(SR_S_function(P48,1),BNP_function(P48))
SR_Test_Function(SR_S_function(P25,1),BNP_function(P25))
SR_Test_Function(SR_S_function(P100,1),BNP_function(P100))

SR_Test_Function(SR_S_function(P10,2),SR_S_function(P10,1))
SR_Test_Function(SR_S_function(P48,2),SR_S_function(P48,1))
SR_Test_Function(SR_S_function(P25,2),SR_S_function(P25,1))
SR_Test_Function(SR_S_function(P100,2),SR_S_function(P100,1))

SR_Test_Function(SR_S_function(P10,2),BNP_function(P10))
SR_Test_Function(SR_S_function(P48,2),BNP_function(P48))
SR_Test_Function(SR_S_function(P25,2),BNP_function(P25))
SR_Test_Function(SR_S_function(P100,2),BNP_function(P100))

SR_Test_Function(SR_S_function(P10,2),PC_MSR_function(P10,k10))
SR_Test_Function(SR_S_function(P48,2),PC_MSR_function(P48,k48))
SR_Test_Function(SR_S_function(P25,2),PC_MSR_function(P25,k25))
SR_Test_Function(SR_S_function(P100,2),PC_MSR_function(P100,k100))

SR_Test_Function(SR_S_function(P10,2),PC_MV_function(P10,k10))
SR_Test_Function(SR_S_function(P48,2),PC_MV_function(P48,k48))
SR_Test_Function(SR_S_function(P25,2),PC_MV_function(P25,k25))
SR_Test_Function(SR_S_function(P100,2),PC_MV_function(P100,k100))

#BN - Simpler ones
SR_Test_Function(BNP_function(P10),R_MV_10)
SR_Test_Function(BNP_function(P48),R_MV_48)
SR_Test_Function(BNP_function(P25),R_MV_25)
SR_Test_Function(BNP_function(P100),R_MV_100)

SR_Test_Function(BNP_function(P10),R_MSR_10)
SR_Test_Function(BNP_function(P48),R_MSR_48)
SR_Test_Function(BNP_function(P25),R_MSR_25)
SR_Test_Function(BNP_function(P100),R_MSR_100)

SR_Test_Function(BNP_function(P10),R_AVP_10)
SR_Test_Function(BNP_function(P48),R_AVP_48)
SR_Test_Function(BNP_function(P25),R_AVP_25)
SR_Test_Function(BNP_function(P100),R_AVP_100)

SR_Test_Function(BNP_function(P10),R_EW_10)
SR_Test_Function(BNP_function(P48),R_EW_48)
SR_Test_Function(BNP_function(P25),R_EW_25)
SR_Test_Function(BNP_function(P100),R_EW_100)


#Severini1 - Simpler ones
SR_Test_Function(SR_S_function(P10,1),R_MV_10)
SR_Test_Function(SR_S_function(P48,1),R_MV_48)
SR_Test_Function(SR_S_function(P25,1),R_MV_25)
SR_Test_Function(SR_S_function(P100,1),R_MV_100)

SR_Test_Function(SR_S_function(P10,1),R_MSR_10)
SR_Test_Function(SR_S_function(P48,1),R_MSR_48)
SR_Test_Function(SR_S_function(P25,1),R_MSR_25)
SR_Test_Function(SR_S_function(P100,1),R_MSR_100)

SR_Test_Function(SR_S_function(P10,1),R_AVP_10)
SR_Test_Function(SR_S_function(P48,1),R_AVP_48)
SR_Test_Function(SR_S_function(P25,1),R_AVP_25)
SR_Test_Function(SR_S_function(P100,1),R_AVP_100)

SR_Test_Function(SR_S_function(P10,1),R_EW_10)
SR_Test_Function(SR_S_function(P48,1),R_EW_48)
SR_Test_Function(SR_S_function(P25,1),R_EW_25)
SR_Test_Function(SR_S_function(P100,1),R_EW_100)

#Severini2 - Simpler ones
SR_Test_Function(SR_S_function(P10,2),R_MV_10)
SR_Test_Function(SR_S_function(P48,2),R_MV_48)
SR_Test_Function(SR_S_function(P25,2),R_MV_25)
SR_Test_Function(SR_S_function(P100,2),R_MV_100)

SR_Test_Function(SR_S_function(P10,2),R_MSR_10)
SR_Test_Function(SR_S_function(P48,2),R_MSR_48)
SR_Test_Function(SR_S_function(P25,2),R_MSR_25)
SR_Test_Function(SR_S_function(P100,2),R_MSR_100)

SR_Test_Function(SR_S_function(P10,2),R_AVP_10)
SR_Test_Function(SR_S_function(P48,2),R_AVP_48)
SR_Test_Function(SR_S_function(P25,2),R_AVP_25)
SR_Test_Function(SR_S_function(P100,2),R_AVP_100)

SR_Test_Function(SR_S_function(P10,2),R_EW_10)
SR_Test_Function(SR_S_function(P48,2),R_EW_48)
SR_Test_Function(SR_S_function(P25,2),R_EW_25)
SR_Test_Function(SR_S_function(P100,2),R_EW_100)

#PCMV - MV
SR_Test_Function(PC_MV_function(P10,k10),R_MV_10)
SR_Test_Function(PC_MV_function(P48,k48),R_MV_48)
SR_Test_Function(PC_MV_function(P25,k25),R_MV_25)
SR_Test_Function(PC_MV_function(P100,k100),R_MV_100)

#PCMSR - MSR
SR_Test_Function(PC_MSR_function(P10,k10),R_MSR_10)
SR_Test_Function(PC_MSR_function(P48,k48),R_MSR_48)
SR_Test_Function(PC_MSR_function(P25,k25),R_MSR_25)
SR_Test_Function(PC_MSR_function(P100,k100),R_MSR_100)

#PCVP - AVP
SR_Test_Function(R_PCVP_10,R_AVP_10)
SR_Test_Function(R_PCVP_48,R_AVP_48)
SR_Test_Function(R_PCVP_25,R_AVP_25)
SR_Test_Function(R_PCVP_100,R_AVP_100)

#LWMSR - MSR
SR_Test_Function(R_LWMSR_10,R_MSR_10)
SR_Test_Function(R_LWMSR_48,R_MSR_48)
SR_Test_Function(R_LWMSR_25,R_MSR_25)
SR_Test_Function(R_LWMSR_100,R_MSR_100)

#LWMV - MV
SR_Test_Function(R_LWMV_10,R_MV_10)
SR_Test_Function(R_LWMV_48,R_MV_48)
SR_Test_Function(R_LWMV_25,R_MV_25)
SR_Test_Function(R_LWMV_100,R_MV_100)

#########################################################################################
########################### Data Description ############################################
#########################################################################################

##Basics
mean(P10)
mean(P48)
mean(P25)
mean(P100)
sd(P10)
sd(P48)
sd(P25)
sd(P100)

##Histograms
par(mfrow=c(2,2))

hist(P10, main="", col="cornflowerblue")
hist(P48, main="", col="cornflowerblue")
hist(P25, main="", col="cornflowerblue")
hist(P100, main="", col="cornflowerblue")

## Skewness - Kurtosis
par(mfrow=c(4,2))

skewP10 = as.matrix(t(skewness(P10)))
KurtoP10 = as.matrix(t(kurtosis(P10)))
mean(skewP10)
mean(KurtoP10)
plot(skewP10, type = 'l', col = 'lightsteelblue', xlab = "Assets", ylab = "Skewness", main = "P10")
plot(KurtoP10, type = 'l', col = 'lightsteelblue',xlab = "Assets", ylab = "Kurtosis", main = "P10")

skewP48 = as.matrix(t(skewness(P48)))
KurtoP48 = as.matrix(t(kurtosis(P48)))
mean(skewP48)
mean(KurtoP48)
plot(skewP48, type = 'l', col = 'lightsteelblue',xlab = "Assets", ylab = "Skewness", main = "P48")
plot(KurtoP48, type = 'l', col = 'lightsteelblue',xlab = "Assets", ylab = "Kurtosis", main = "P48")

skewP25 = as.matrix(t(skewness(P25)))
KurtoP25 = as.matrix(t(kurtosis(P25)))
mean(skewP25)
mean(KurtoP25)
plot(skewP25, type = 'l', col = 'lightsteelblue',xlab = "Assets", ylab = "Skewness", main = "P25")
plot(KurtoP25, type = 'l', col = 'lightsteelblue',xlab = "Assets", ylab = "Kurtosis", main = "P25")

skewP100 = as.matrix(t(skewness(P100)))
KurtoP100 = as.matrix(t(kurtosis(P100)))
mean(skewP100)
mean(KurtoP100)
plot(skewP100, type = 'l', col = 'lightsteelblue',xlab = "Assets", ylab = "Skewness", main = "P100")
plot(KurtoP100, type = 'l', col = 'lightsteelblue',xlab = "Assets", ylab = "Kurtosis", main = "P100")


# Visualize correlation matrix
corrplot(Qtest, method = "color", type = "upper", diag = FALSE, tl.col = "black", addCoef.col = "black")

# Hierarchical clustering
hc <- hclust(dist(cor(P10)))
dend <- as.dendrogram(hc)

# Visualize dendrogram
plot(dend, main = "Dendrogram of Asset Correlations", ylab = "Euclidean Distance", xlab = "Assets")

# Create the dataset
returns_density <- A$density
data <- data.frame(return_density = returns_density)

ggplot(data, aes(x = return_density)) +
  stat_ecdf(geom = "step", color = "steelblue") +
  labs(title = "Cumulative Distribution Function (CDF) Plot of Return Densities", x = "Return Density", y = "Cumulative Probability") +
  theme_minimal()


