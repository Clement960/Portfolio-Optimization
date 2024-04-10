#########################################################################################
######################### Equally-weighted portfolio ####################################
#########################################################################################
#------------
#Date: 09/04/2024
#
#
#------------

  equally_weighted = function(data){
    
    #Weights
    if (10 <= nrow(data)/6 && nrow(data)/6  < 20){
      estperiods = 10
      periods = nrow(data)/6 - estperiods
      } else if (nrow(data)/6 >= 20){
        estperiods = 20
        periods = nrow(data)/6 - estperiods
        } else {stop("Dataset is to small, must be at least 60 months.")}
    
    W = matrix(1/ncol(data),  ncol(data))
    
    
    #Returns
    R = array(dim=c(6,periods))
    
    a=0
    b=0
    for (i in 1:periods){
      a[i]=estperiods*6 -5+i*6
      b[i]=estperiods*6+i*6
      
      R[,i] =(data[a[i]:b[i],]%*%W)
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
