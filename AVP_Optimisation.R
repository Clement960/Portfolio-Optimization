

###OPTIMISATION METHOD - Asset-Variance-Parity Portfolio
w_function = function(data){
  cov_matrix = cov(data)
  
  # Define the objective function
  obj_func <- function(w, cov_matrix) {
    return(w %*% cov_matrix %*% w - sum(log(w)))
  }
  
  # Define the initial value for w
  n_assets <- ncol(cov(data))
  w_init <- rep(1/n_assets, n_assets)
  
  # Set the optimization control parameters
  ctrl <- list(fnscale = 1,  # minimise the objective function
               factr = 1e-8,  # absolute tolerance for convergence
               maxit = 1000)  # maximum number of iterations
  
  # Optimize the objective function subject to the constraint
  result <- optim(w_init, obj_func, cov_matrix = cov_matrix, method = "L-BFGS-B", lower = rep(0, n_assets), control = ctrl)
  
  # Extract the optimal solution
  w_optimal <- result$par/sum(result$par)
  return(w_optimal)
}
