#### Import all the codes ####
#setwd("C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/github/master_thesis")
source("source_file_common_all.R")
source("script_CVM.R")
source("script_CVMSE.R")
source("script_lepskichen.R")
source("script_lepskiboot.R")



#### Fonction ####
J_opt_data_fixed <- function(simul_data, J_CV_bs, J_CV_ns, p_train, x_evaluation, degree, n_boot_lepski){
  W <- simul_data$W
  Y <- simul_data$Y
  Z <- simul_data$Z 
  
  m_m = max(Z) - min(Z)
    
  length_J_CV <- length(J_CV_bs)
  
  #compute all gamma_J on the data 
  list_M_boot_bs <- vector("list", length_J_CV)
  list_M_boot_ns <- vector("list", length_J_CV)
  list_gamma_bs <- vector("list", length_J_CV)
  list_gamma_ns <- vector("list", length_J_CV)
  
  for (j_index in 1:length_J_CV){
    J_bs <- J_CV_bs[j_index]
    J_ns <- J_CV_ns[j_index]
    M_boot_j_bs <- compute_M_bootstrap(J_bs, W, Z, Y, degree, create_dyadic_P_splines_bs) 
    list_M_boot_bs[[j_index]] <- M_boot_j_bs
    list_gamma_bs[[j_index]] = M_boot_j_bs%*%Y
    
    M_boot_j_ns <- compute_M_bootstrap(J_ns, W, Z, Y, degree, create_dyadic_P_splines_ns)
    list_M_boot_ns[[j_index]] <- M_boot_j_ns
    list_gamma_ns[[j_index]] = M_boot_j_ns%*%Y
  }
  
  #compute the g_hat_J for all J on the data on the values of x_evaluation
  list_g_hat_J_bs <- vector("list", length_J_CV)
  list_g_hat_J_ns <- vector("list", length_J_CV)
  for (j_index in 1:length_J_CV){
    J_bs <- J_CV_bs[j_index]
    J_ns <- J_CV_ns[j_index]
    basis_b_spline <- create_dyadic_P_splines_bs(x_evaluation, Z, J_bs, degree)
    basis_n_spline <- create_dyadic_P_splines_ns(x_evaluation, Z, J_ns, degree)
    list_g_hat_J_bs[[j_index]] <- basis_b_spline %*% list_gamma_bs[[j_index]]
    list_g_hat_J_ns[[j_index]] <- basis_n_spline %*% list_gamma_ns[[j_index]]
    #on a calculé les g_hat_J pour tous les J -> on a besoin de renvoyer que les indices et cette liste
  }
  
  #cross validation on M and MSE
  zero_indices_bs = which(sapply(list_g_hat_J_bs, function(x) all(x == 0))) #filter indices
  zero_indices_ns = which(sapply(list_g_hat_J_ns, function(x) all(x == 0)))
  
  if (length(zero_indices_bs)>0){
    new_J_CV_bs <- J_CV_bs[-zero_indices_bs]}
  else{
    new_J_CV_bs <- J_CV_bs}
  if (length(zero_indices_ns)>0){
    new_J_CV_ns <- J_CV_ns[-zero_indices_ns]}
  else{
    new_J_CV_ns <- J_CV_ns
  }
  
  #do the splitting for cross validation : to have the same data for both algorithms
  sampled_data <- sample_train_test(Z, Y, W, p_train)
  Z_train <- sampled_data$Z_train
  Y_train <- sampled_data$Y_train
  W_train <- sampled_data$W_train
  Z_validation <- sampled_data$Z_validation
  Y_validation <- sampled_data$Y_validation
  W_validation <- sampled_data$W_validation
  
  
  if (length(new_J_CV_bs)>1){
    J_opt_CV_M_bs <- optimization_CV_M(Z_train, W_train, Y_train, Z_validation, W_validation, Y_validation,
                                       new_J_CV_bs, p_train, 
                                       degree, x_evaluation, create_dyadic_P_splines_bs)
    
    J_opt_CVMSE_bs <- optimization_CV_MSE(Z_train, W_train, Y_train, Z_validation, W_validation, Y_validation,
                                          new_J_CV_bs, p_train, 
                                          degree, x_evaluation, create_dyadic_P_splines_bs)
  }
  else if(length(new_J_CV_bs)==1){ #si on a qu'un seul choix
    J_opt_CV_M_bs <- new_J_CV_bs
    J_opt_CVMSE_bs <- new_J_CV_bs
    }
  else{
    J_opt_CV_M_bs <- 0
    J_opt_CVMSE_bs <- 0
  }
  
  if (length(new_J_CV_ns)>1){
    J_opt_CV_M_ns <- optimization_CV_M(Z_train, W_train, Y_train, Z_validation, W_validation, Y_validation,
                                       new_J_CV_ns, p_train, 
                                       degree, x_evaluation, create_dyadic_P_splines_ns)
    
    J_opt_CVMSE_ns <- optimization_CV_MSE(Z_train, W_train, Y_train, Z_validation, W_validation, Y_validation,
                                          new_J_CV_ns, p_train, 
                                          degree, x_evaluation, create_dyadic_P_splines_ns)
  }
  else if (length(new_J_CV_ns)==1){
    J_opt_CV_M_ns <- new_J_CV_ns
    J_opt_CVMSE_ns <- new_J_CV_ns
  }
  else{
    J_opt_CV_M_ns <- 0
    J_opt_CVMSE_ns <- 0
  }
  
  rm(sampled_data, Z_train, Z_validation, Y_train, Y_validation, W_train, W_validation)
  gc()
  
  #pour lepski
  bs_bool = 1
  J_max_bs = compute_J_max(W, Z, Y, degree, bs_bool)
  if(J_max_bs > max(J_CV_bs)){
    J_max_bs = max(J_CV_bs)
  }
  I_hat_bs = seq(as.integer(0.1 * (log(J_max_bs)^2)), as.integer(J_max_bs), by = 1) #on commence au moins à 2 parce que 1 c'est pas une bonne solution
  I_hat_bs = sort(I_hat_bs[sapply(I_hat_bs,valid_dim_b_splines, degree)]) #select only the valid dimensions
  
  bs_bool = 0
  J_max_ns = compute_J_max(W, Z, Y, degree, bs_bool)
  if(J_max_ns > max(J_CV_ns)){
    J_max_ns = max(J_CV_ns)
  }
  I_hat_ns = seq(as.integer(0.1 * (log(J_max_ns)^2)), as.integer(J_max_ns), by = 1) #on commence au moins à 2 parce que 1 c'est pas une bonne solution
  I_hat_ns = sort(I_hat_ns[sapply(I_hat_ns,valid_dim_ns, degree)]) #select only the valid dimensions
  
  # filtrer les valeurs possibles : on enlève ceux où gamma est zéro
  if (length(zero_indices_bs)>0){
    I_hat_bs_new <- I_hat_bs[-zero_indices_bs]
  }
  else{
    I_hat_bs_new <- I_hat_bs
  }
  if (length(zero_indices_ns)>0){
    I_hat_ns_new <- I_hat_ns[-zero_indices_ns]
  }
  else{
    I_hat_ns_new <- I_hat_ns
  }
  

  if (length(I_hat_bs_new) > 1){
    c_0 <- 10
    J_opt_lepski_bs <- lepski_chen(I_hat_bs_new, c_0, m_m, list_g_hat_J_bs, W, Z, Y, degree, 
                                   valid_dim_b_splines, create_dyadic_P_splines_bs)
    
    # lepski boot
    J_opt_lepskiboot_bs <- lepski_bootstrap(I_hat_bs_new, n_boot_lepski, list_gamma_bs, list_M_boot_bs, 
                                            list_g_hat_J_bs, valid_dim_b_splines, x_evaluation, 
                                            W, Z, Y, degree, create_dyadic_P_splines_bs)
    

    }
  else if (length(I_hat_bs_new)==1){
    J_opt_lepski_bs <- I_hat_ns_new
    J_opt_lepskiboot_bs <- I_hat_ns_new
  }
  else{
    J_opt_lepski_bs <- 0
    J_opt_lepskiboot_bs <- 0
  }
  
  if (length(I_hat_ns_new)>1){
    J_opt_lepski_ns <- lepski_chen(I_hat_ns_new, c_0, m_m, list_g_hat_J_ns, W, Z, Y, degree, 
                                   valid_dim_ns, create_dyadic_P_splines_ns)

    J_opt_lepskiboot_ns <- lepski_bootstrap(I_hat_ns_new, n_boot_lepski, list_gamma_ns, list_M_boot_ns,
                                            list_g_hat_J_ns, valid_dim_ns, x_evaluation, 
                                            W, Z, Y, degree, create_dyadic_P_splines_ns)
    }
  else if (length(I_hat_ns_new)==1){
    J_opt_lepski_ns <- I_hat_ns_new
    J_opt_lepskiboot_ns <- I_hat_ns_new
  }
  else{
    J_opt_lepski_ns <- 0
    J_opt_lepskiboot_ns <- 0
  }
  
  rm(list_M_boot_bs, list_M_boot_ns, I_hat_bs, I_hat_ns)
  gc()
  
  
  #return the results
  return(list(J_opt_CV_M_bs = J_opt_CV_M_bs, J_opt_CV_M_ns = J_opt_CV_M_ns,
              J_opt_CVMSE_bs = J_opt_CVMSE_bs, J_opt_CVMSE_ns = J_opt_CVMSE_ns, 
              J_opt_lepski_bs = J_opt_lepski_bs, J_opt_lepski_ns = J_opt_lepski_ns,
              J_opt_lepskiboot_bs = J_opt_lepskiboot_bs, J_opt_lepskiboot_ns = J_opt_lepskiboot_ns,
              list_gamma_bs = list_gamma_bs, list_gamma_ns = list_gamma_ns,
              list_g_hat_J_bs = list_g_hat_J_bs, list_g_hat_J_ns = list_g_hat_J_ns,
              W = W , Y = Y, Z = Z))
  }


# 
# J_bs_ <- c(4, 5, 7, 11, 19)
# J_ns_ <- c(2, 3, 5, 9, 17)
# p = 0.5
# x_eval = seq(-2, 2, length.out = 100)
# deg = 3
# n_boot = 10
# simul <- simulate_data_3(c(1000, 0.5, 0.9), g_sim_3, 2)
# # 
# # #setwd("C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/simulation_results/4_final")
# # load("data_2000_rhozw0.9_rhouv0.8_case3_n1000.R")
# # 
#  test <- J_opt_data_fixed(simul, J_bs_, J_ns_, p, x_eval, deg, n_boot)
# 
# 
# simul_all[[2]]$W
# simul$W









