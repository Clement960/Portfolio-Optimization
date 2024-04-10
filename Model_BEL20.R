#Library-------------------------------------------------------------------
pacman::p_load(netstat, RSelenium, rio,rvest,magrittr,dtplyr,forecast,DT,zoo,lubridate,hrbrthemes, data.table,arrow, seasonal,seasonalview,readxl, tools, xts,glmnet,dplyr,POET,stcov,fastICA,nlshrink,rmgarch,moments, PerformanceAnalytics, quadprog,NMOF, riskParityPortfolio, tidyquant,timetk,ggplot2,foreach,MASS,CovTools,readr,paran,factoextra,tidyverse,scales,tidyr,corrplot,reshape2,gridExtra,palmerpenguins,SciViews,dfms,boot,PeerPerformance,ggthemes)

source("J:/Minimum_Variance.R")
source("J:/Mean_Variance.R")
source("J:/Equally_Weighted.R")
source("J:/Asset_Variance_Parity.R")
source("J:/Ledoit-Wolf_MV-MSR.R")

#Data Importation--------------------------------------------------------------------
setwd("J:/C_Piedboeuf/Doc Admin/ING/Business Panel/Data")



# # URL of the Investing.com page
# url <- "https://fr.investing.com/indices/bel-20-historical-data"
# 
# # Read the HTML content of the webpage
# webpage <- read_html(url)
# 
# # Extract the historical data table using CSS selector
# historical_data <- webpage %>%
#   html_node("#curr_table") %>%
#   html_table()



#Data source:--------------------------------------------------------------------
#https://fr.investing.com/indices/bel-20-historical-data

BEL20 <- read_csv("BEL 20 - Données Historiques (1).csv")[11:167,3] 
ABI   <- read_csv("ABI - Données Historiques.csv")[1:97,3]
ACKB  <- read_csv("ACKB - Données Historiques.csv")[1:97,3] 
AGES  <- read_csv("AGES - Données Historiques.csv")[1:97,3] 
AOO   <- read_csv("AOO - Données Historiques.csv")[1:97,3]  
APAM  <- read_csv("APAM - Données Historiques.csv")[1:97,3] 
ARGX  <- read_csv("ARGX - Données Historiques.csv")[1:97,3] 
BAR   <- read_csv("BAR - Données Historiques.csv")[1:97,3]  
COFB  <- read_csv("COFB - Données Historiques.csv")[1:97,3] 
ELI   <- read_csv("ELI - Données Historiques.csv")[1:97,3]  
GBLB  <- read_csv("GBLB - Données Historiques.csv")[1:97,3] 
GLPG  <- read_csv("GLPG - Données Historiques.csv")[1:97,3] 
IETB  <- read_csv("IETB - Données Historiques.csv")[1:97,3] 
KBC   <- read_csv("KBC - Données Historiques.csv")[1:97,3]  
MLXS  <- read_csv("MLXS - Données Historiques.csv")[1:97,3] 
PROX  <- read_csv("PROX - Données Historiques.csv")[1:97,3] 
SOF   <- read_csv("SOF - Données Historiques.csv")[1:97,3]  
SOLB  <- read_csv("SOLB - Données Historiques.csv")[1:97,3] 
UCB   <- read_csv("UCB - Données Historiques.csv")[1:97,3]  
UMI   <- read_csv("UMI - Données Historiques.csv")[1:97,3]  
WDPP  <- read_csv("WDPP - Données Historiques.csv")[1:97,3] 

Synthetic_BEL20 <- cbind(ABI,ACKB, AGES, AOO, APAM, ARGX, BAR, COFB, ELI, GBLB, GLPG, IETB, KBC, MLXS, PROX, SOF, SOLB, UCB, UMI, WDPP )
column_names <- c("ABI", "ACKB", "AGES", "AOO", "APAM", "ARGX", "BAR", "COFB", "ELI", "GBLB", "GLPG", "IETB", "KBC", "MLXS", "PROX", "SOF", "SOLB", "UCB", "UMI", "WDPP")
names(Synthetic_BEL20) <- column_names


# Calculate the percentage change for each column
Synthetic_BEL20 <- (apply(Synthetic_BEL20, 2, function(x) {
  diff(x) / x[-length(x)]
}))



#Minimum-variance Portfolio:--------------------------------------------------------------------
out_MV = minimum_variance(Synthetic_BEL20)

#Mean-variance Portfolio:--------------------------------------------------------------------
out_MSR = mean_variance(Synthetic_BEL20)

#Equally-Weighted Portfolio:--------------------------------------------------------------------
out_EW = equally_weighted(Synthetic_BEL20)

#Asset-Variance-Parity Portfolio:--------------------------------------------------------------------
out_AVP = asset_variance_parity(Synthetic_BEL20)

#Ledoit-Wolf Minimum-Variance Portfolio:--------------------------------------------------------------------
out_AVP = Ledoit_wolf_MV(Synthetic_BEL20)

#Ledoit-Wolf Mean-Variance Portfolio:--------------------------------------------------------------------
out_AVP = Ledoit_wolf_MSR(Synthetic_BEL20)





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


