#########################################################################################
######################### Asset-Variance-Parity portfolio ###############################
#########################################################################################

source("J:/AVP_Optimisation.R")


asset_variance_parity = function(data){
  
  # #PACKAGE METHOD
  # W= riskParityPortfolio(cov(data))$w     
  
  
  #Risk proportion of every asset
  RP =array(dim = c(1,ncol(data)))
  W_optimal = w_function(data)
  for (i in 1:ncol(data)){
    RP[i] =(W_optimal[i] %*% (cov(data) %*% W_optimal)[i])/(W_optimal %*% (cov(data) %*% W_optimal))
    }
  
  #Weights (Computed in sample)
  if (10 <= nrow(data)/6 && nrow(data)/6  < 20){
    estperiods = 10
    periods = nrow(data)/6 - estperiods
    W = array(dim=c(ncol(data),periods))
    } else if (nrow(data)/6 >= 20){
      estperiods = 20
      periods = nrow(data)/6 - estperiods
      W = array(dim=c(ncol(data),periods))
      } else {stop("Dataset is to small, must be at least 60 months.")}
  
  
  a=0
  b=0
  for(i in 1:periods){
    a[i]=i*6-5
    b[i]= estperiods*6 -5 +i*6
    
    W[,i] = t(w_function(data[a[i]:b[i],]))
    }
  
  #Returns
  R = array(dim=c(6,periods))
  
  
  for (i in 1:periods){
    a[i]=estperiods*6 -5+i*6
    b[i]=estperiods*6+i*6
    
    R[,i]<-(data[a[i]:b[i],]%*%W[,i])
    }
  
  out <- list(
    AR = mean(R)*12,
    SD = sd(R)*sqrt(12),
    SR = (mean(R)*12)/(sd(R)*sqrt(12)),
    W = W,
    R = R
  )
  
  return(out)
}
