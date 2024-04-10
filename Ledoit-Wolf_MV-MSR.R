#########################################################################################
######################### Ledoit-Wolf Covariance Estimator ##############################
#########################################################################################

# LW mean-variance:--------------------------------------------------------------------

Ledoit_wolf_MSR = function(data){
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
  
    "W_LWMSR_10[,i] = (solve(CovEst.2003LW((P10[a[i]:b[i],]))$S)%*%colMeans(P10[a[i]:b[i],]))/as.numeric(t(array(1,dim=c(10,1)))%*%solve(CovEst.2003LW((P10[a[i]:b[i],]))$S)%*%colMeans(P10[a[i]:b[i],]))"
  
    W[,i] = (solve(linshrink_cov((data[a[i]:b[i],])))%*%colMeans(data[a[i]:b[i],]))/as.numeric(t(array(1,dim=c(ncol(data),1)))%*%solve(linshrink_cov((data[a[i]:b[i],])))%*%colMeans(data[a[i]:b[i],]))
  
  }

  #Returns (Out of sample)
  R = array(dim=c(6,periods))
  
  for (i in 1:periods){
    a[i]=estperiods*6 -5+i*6
    b[i]=estperiods*6+i*6
    
    R[,i] = (data[a[i]:b[i],]%*%W[,i])
  }
  
  out <- list(
    AR = mean(R)*12,
    SD = sd(R)*sqrt(12),
    SR = (mean(R)*12)/(sd(R)*sqrt(12)),
    W = W,
    R = R
    
  )

boxplot = boxplot(t(W), main = "Ledoit-Wolf Mean-Variance")


}

# LW mean-variance:--------------------------------------------------------------------

Ledoit_wolf_MV = function(data){
  
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
    
    W[,i] = minvar(linshrink_cov(((data[a[i]:b[i],]))))
    
  }
  
  #Returns (Out of sample)
  R = array(dim=c(6,periods))
  
  for (i in 1:periods){
    a[i]=estperiods*6 -5+i*6
    b[i]=estperiods*6+i*6
    
    R[,i] = (data[a[i]:b[i],]%*%W[,i])
    }
  
  out <- list(
    AR = mean(R)*12,
    SD = sd(R)*sqrt(12),
    SR = (mean(R)*12)/(sd(R)*sqrt(12)),
    W = W,
    R = R
  )

    boxplot(t(W), main = "Ledoit-Wolf Mean-Variance",ylim = c(0, 0.5), col = "pink1")
}

