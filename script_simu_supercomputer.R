#### Import all the codes ####
setwd("C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/github/master_thesis")
source("source_file_common_all.R")
source("script_CVM.R")
source("script_CVMSE.R")
source("script_lepskichen.R")
#source("script_lepskiboot.R")

#### Fonction ####
J_opt_data_fixed <- function(simul_data, J_CV_bs, J_CV_ns, p_train, x_evaluation, degree){
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
    list_gamma_ns[[j_index]] = M_boot_j_ns%*%Y}
  
  
  #compute the g_hat_J for all J on the data on the values of x_evaluation
  list_g_hat_J_bs <- vector("list", length_J_CV)
  list_g_hat_J_ns <- vector("list", length_J_CV)
  for (j_index in 1:length_J_CV){
    J_bs <- J_CV_bs[j_index]
    J_ns <- J_CV_ns[j_index]
    basis_b_spline <- create_dyadic_P_splines_bs(x_evaluation, Z, J_bs, degree)
    basis_n_spline <- create_dyadic_P_splines_ns(x_evaluation, Z, J_ns, degree)
    list_g_hat_J_bs[[j_index]] <- basis_b_spline %*% list_gamma_bs[[j_index]]
    list_g_hat_J_ns[[j_index]] <- basis_n_spline %*% list_gamma_ns[[j_index]]}
  
  #filter for cross-validation the indices possible
  zero_indices_bs = which(sapply(list_g_hat_J_bs, function(x) all(x == 0)))
  zero_indices_ns = which(sapply(list_g_hat_J_ns, function(x) all(x == 0)))
  
  new_J_CV_bs <- J_CV_bs[-zero_indices_bs]
  new_J_CV_ns <- J_CV_ns[-zero_indices_ns]
  
  #cross validation on M
  J_opt_CV_M_bs <- optimization_CV_M(Z, W, Y, new_J_CV_bs, p_train, 
                                     degree, x_evaluation, create_dyadic_P_splines_bs)
  j_index_ <- which(J_CV_bs == J_opt_CV_M_bs)
  g_hat_J_bs_CVM = list_g_hat_J_bs[[j_index_]]
  
  J_opt_CV_M_ns <- optimization_CV_M(Z, W, Y, new_J_CV_ns, p_train, 
                                     degree, x_evaluation, create_dyadic_P_splines_ns)
  j_index_ <- which(J_CV_ns == J_opt_CV_M_ns)
  g_hat_J_ns_CVM = list_g_hat_J_ns[[j_index_]]
  
  #cross validation on MSE
  J_opt_CVMSE_bs <- optimization_CV_MSE(Z, W, Y, new_J_CV_bs, p_train, 
                                       degree, x_evaluation, create_dyadic_P_splines_bs)
  j_index_ <- which(J_CV_bs == J_opt_CVMSE_bs)
  g_hat_J_bs_CVMSE = list_g_hat_J_bs[[j_index_]]
  
  J_opt_CVMSE_ns <- optimization_CV_MSE(Z, W, Y, new_J_CV_ns, p_train, 
                                       degree, x_evaluation, create_dyadic_P_splines_ns)
  j_index_ <- which(J_CV_ns == J_opt_CVMSE_ns)
  g_hat_J_ns_CVMSE = list_g_hat_J_ns[[j_index_]]
  
  #pour lepski
  bs_bool = 1
  J_max_bs = compute_J_max(W, Z, Y, degree, bs_bool)
  I_hat_bs = seq(as.integer(0.1 * (log(J_max_bs)^2)), as.integer(J_max_bs), by = 1) #on commence au moins à 2 parce que 1 c'est pas une bonne solution
  I_hat_bs = sort(I_hat_bs[sapply(I_hat_bs,valid_dim_b_splines, degree)]) #select only the valid dimensions
  
  bs_bool = 0
  J_max_ns = compute_J_max(W, Z, Y, degree, bs_bool)
  I_hat_ns = seq(as.integer(0.1 * (log(J_max_ns)^2)), as.integer(J_max_ns), by = 1) #on commence au moins à 2 parce que 1 c'est pas une bonne solution
  I_hat_ns = sort(I_hat_ns[sapply(I_hat_ns,valid_dim_ns, degree)]) #select only the valid dimensions
  
  #lepski 1 : attention mettre une condition pour pas sélectionner un gamma = 0 
  c_0 <- 10
  J_opt_lepski_bs <- lepski_chen(I_hat_bs, c_0, m_m, list_g_hat_J_bs, W, Z, Y, degree, 
                                 valid_dim_b_splines, create_dyadic_P_splines_bs)
  j_index_ <- which(J_CV_bs == J_opt_lepski_bs)
  g_hat_J_bs_lepski = list_g_hat_J_bs[[j_index_]]
  
  J_opt_lepski_ns <- lepski_chen(I_hat_ns, c_0, m_m, list_g_hat_J_ns, W, Z, Y, degree, 
                                 valid_dim_ns, create_dyadic_P_splines_ns)
  j_index_ <- which(J_CV_ns == J_opt_lepski_ns)
  g_hat_J_ns_lepski = list_g_hat_J_ns[[j_index_]]
  
  #lepski boot
  #J_opt_lepskiboot_bs <- lepski_bootstrap(n_boot, valid_dim_b_splines, x_evaluation, 
                                          #W, Z, Y, degree, create_dyadic_P_splines_bs)
  #j_index_ <- which(J_CV_bs == J_opt_lepskiboot_bs)
  #g_hat_J_bs_lepski_boot = list_g_hat_J_bs[[j_index_]]
  
  #J_opt_lepskiboot_ns <- lepski_bootstrap(n_boot, valid_dim_ns, x_evaluation, 
                                          #W, Z, Y, degree, create_dyadic_P_splines_ns)
  #j_index_ <- which(J_CV_bs == J_opt_lepskiboot_ns)
  #g_hat_J_ns_lepski_boot = list_g_hat_J_ns[[j_index_]]
  
  
  #return the results
  return(list(J_opt_CV_M_bs = J_opt_CV_M_bs, J_opt_CV_M_ns = J_opt_CV_M_ns,
              J_opt_CVMSE_bs = J_opt_CVMSE_bs, J_opt_CVMSE_ns = J_opt_CVMSE_ns, 
              J_opt_lepski_bs = J_opt_lepski_bs, J_opt_lepski_ns = J_opt_lepski_ns,
              #J_opt_lepskiboot_bs = J_opt_lepskiboot_bs, J_opt_lepskiboot_ns = J_opt_lepskiboot_ns,
              list_gamma_bs = list_gamma_bs, list_gamma_ns = list_gamma_ns,
              g_hat_J_bs_CVM = g_hat_J_bs_CVM, g_hat_J_ns_CVM = g_hat_J_ns_CVM,
              g_hat_J_bs_CVMSE = g_hat_J_bs_CVMSE, g_hat_J_ns_CVMSE = g_hat_J_ns_CVMSE,
              g_hat_J_bs_lepski = g_hat_J_bs_lepski, g_hat_J_ns_lepski = g_hat_J_ns_lepski#,
              #g_hat_J_bs_lepski_boot = g_hat_J_bs_lepski_boot, g_hat_J_ns_lepski_boot = g_hat_J_ns_lepski_boot
              ))
  
}

#simul <- simulate_data_3(c(200, 0.5, 0.9), g_sim_3, 2)
#J_bs_ <- c(5, 7, 11, 19, 35)
#J_ns_ <- c(3, 5, 9, 17, 33)
#p = 0.5
#x_eval = seq(-2, 2, length.out = 100)
#deg = 3


#test <- J_opt_data_fixed(simul, J_bs_, J_ns_, p, x_eval, deg)



#### MC procedure TO DO ####

list_W <- vector("list", n_MC)
list_Z <- vector("list", n_MC)
list_Y <- vector("list", n_MC)


for (n in 1:n_MC){
  #generate data 
  simul <- simulate_data_3(data_param, g_0, case)
  W <- simul$W
  Y <- simul$Y
  Z <- simul$Z
  
  list_W[[n]] <- W
  list_Y[[n]] <- Y
  list_Z[[n]] <- Z
  
  
  
  }

