#### Selection of J by cross-validation on the MSE ####
calcul_M_g_hat <- function(gamma_hat, J, Z, Y, W, degree, create_P){ #ok testé
  Omega <- create_W(W)
  a = min(Z)
  b = max(Z)
  n = length(Z)
  g_hat_on_Z <- rep(0, n)
  basis <- create_P(Z, Z, J, degree)
  g_hat_on_Z <- basis%*%gamma_hat
  for (i in 1:n){
    for (j in 1:n){
      M = M + (Y[i] - g_hat_on_Z[i])*(Y[j] - g_hat_on_Z[j])*Omega[i,j]
    }
  }
  rm(Omega, g_hat_on_Z, d, basis)
  gc()
  return(M/(n**2))}

calcul_M_true_g <- function(g, J, Z, Y, W, case){
  M = 0
  n = length(Z)
  Omega <- create_W(W)
  g_hat <- sapply(Z, g, case = case)
  
  for (i in 1:n){
    for (j in 1:n){
      M = M + (Y[i] - g_hat[i])*(Y[j] - g_hat[j])*Omega[i,j]
    }}
  
  rm(Omega, d, M, g_hat)
  gc()
  return(M/(n**2))}


sample_train_test <- function(Z, Y, W, p_train){
  df <- data.frame(Z = Z, Y = Y, W = W)
  
  # Create bins based on the X variable
  bins <- cut(df$Z, breaks = 3, labels = FALSE)
  index_train <- createDataPartition(bins, p = p_train, list = FALSE)
  index_val <- setdiff(seq_len(nrow(df)), index_train)
  
  train_set <- df[index_train, ]
  val_set   <- df[index_val, ]
  
  return(list(
    Z_train = train_set$Z, Y_train = train_set$Y, W_train = train_set$W,
    Z_validation = val_set$Z, Y_validation = val_set$Y, W_validation = val_set$W
  ))
}

optimization_CV_MSE <- function(Z_train, W_train, Y_train, Z_validation, W_validation, Y_validation,
                                vect_J_to_test, p_train, degree, x_grid, create_P){
  gamma_hat_list <- vector("list", length(vect_J_to_test))
  #a = min(Z)
  #b = max(Z)
  MSE_values <- numeric(length(vect_J_to_test))
  #n = length(Z)
  
  n_train <- length(Z_train)
  
  # compute gamma_hat_train and g_hat_train for the training set
  for (j_index in 1:length(vect_J_to_test)){
    J = vect_J_to_test[j_index]
    #compute gamma_hat_train
    gamma_hat <- estimation_gamma(J,W_train,Z_train,Y_train, degree, create_P)
    gamma_hat_list[[j_index]] <- gamma_hat
    
    # compute the prediction of g_hat_train on the test set 
    n_val <- length(Z_validation)
    basis <- create_P(Z_validation, Z_train, J, degree)
    g_hat_on_Z_val <-  basis%*%gamma_hat
    
    # compute MSE and store it
    MSE = sum((g_hat_on_Z_val - Z_validation)**2)
    MSE_values[j_index] = MSE
  }
  
  # return J corresponding to the smallest MSE
  J_opt = vect_J_to_test[which.min(MSE_values)]
  
  rm(Z_train, Z_validation, Y_train, Y_validation, W_train, W_validation, 
     gamma_hat_list, MSE_values)
  gc()
  
  #return(list(J_opt = J_opt, gamma_hat_j_opt = gamma_hat, g_hat_on_x_opt = g_hat_on_x_opt, gamma_hat_list = gamma_hat_list, MSE_values = MSE_values))
  return(J_opt)
}

#done