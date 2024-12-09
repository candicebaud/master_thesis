#### Import all the codes ####
setwd("C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/github/master_thesis")
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
    list_g_hat_J_ns[[j_index]] <- basis_n_spline %*% list_gamma_ns[[j_index]]
    #on a calculé les g_hat_J pour tous les J -> on a besoin de renvoyer que les indices et cette liste
  }
  
  
  #cross validation on M and MSE
  zero_indices_bs = which(sapply(list_g_hat_J_bs, function(x) all(x == 0))) #filter indices
  zero_indices_ns = which(sapply(list_g_hat_J_ns, function(x) all(x == 0)))
  
  new_J_CV_bs <- J_CV_bs[-zero_indices_bs] #nouvelle liste d'indices
  new_J_CV_ns <- J_CV_ns[-zero_indices_ns]
  
  #do the splitting for cross validation : to have the same data for both algorithms
  sampled_data <- sample_train_test(Z, Y, W, p_train)
  Z_train <- sampled_data$Z_train
  Y_train <- sampled_data$Y_train
  W_train <- sampled_data$W_train
  Z_validation <- sampled_data$Z_validation
  Y_validation <- sampled_data$Y_validation
  W_validation <- sampled_data$W_validation
  
  if (length(new_J_CV_bs)>0){
    J_opt_CV_M_bs <- optimization_CV_M(Z_train, W_train, Y_train, Z_validation, W_validation, Y_validation,
                                       new_J_CV_bs, p_train, 
                                       degree, x_evaluation, create_dyadic_P_splines_bs)
    #j_index_ <- which(J_CV_bs == J_opt_CV_M_bs)
    #g_hat_J_bs_CVM = list_g_hat_J_bs[[j_index_]]
    
    J_opt_CVMSE_bs <- optimization_CV_MSE(Z_train, W_train, Y_train, Z_validation, W_validation, Y_validation,
                                          new_J_CV_bs, p_train, 
                                          degree, x_evaluation, create_dyadic_P_splines_bs)
    #j_index_ <- which(J_CV_bs == J_opt_CVMSE_bs)
    #g_hat_J_bs_CVMSE = list_g_hat_J_bs[[j_index_]]
    }
  else{
    J_opt_CV_M_bs <- 0
    #g_hat_J_bs_CVM <- rep(0, length(x_evaluation))
    J_opt_CVMSE_bs <- 0
    #g_hat_J_bs_CVMSE <- rep(0, length(x_evaluation))
  }
  
  
  if (length(new_J_CV_ns)>0){
    J_opt_CV_M_ns <- optimization_CV_M(Z_train, W_train, Y_train, Z_validation, W_validation, Y_validation,
                                       new_J_CV_ns, p_train, 
                                       degree, x_evaluation, create_dyadic_P_splines_ns)
    #j_index_ <- which(J_CV_ns == J_opt_CV_M_ns)
    #g_hat_J_ns_CVM = list_g_hat_J_ns[[j_index_]]
    
    J_opt_CVMSE_ns <- optimization_CV_MSE(Z_train, W_train, Y_train, Z_validation, W_validation, Y_validation,
                                          new_J_CV_ns, p_train, 
                                          degree, x_evaluation, create_dyadic_P_splines_ns)
    #j_index_ <- which(J_CV_ns == J_opt_CVMSE_ns)
    #g_hat_J_ns_CVMSE = list_g_hat_J_ns[[j_index_]]}
  }
  else{
    J_opt_CV_M_ns <- 0
    #g_hat_J_ns_CVM <- rep(0, length(x_evaluation))
    J_opt_CVMSE_ns <- 0
    #g_hat_J_ns_CVMSE <- rep(0, length(x_evaluation))
  }
  
  rm(sampled_data, Z_train, Z_validation, Y_train, Y_validation, W_train, W_validation)
  gc()
  
  #pour lepski
  bs_bool = 1
  J_max_bs = compute_J_max(W, Z, Y, degree, bs_bool)
  I_hat_bs = seq(as.integer(0.1 * (log(J_max_bs)^2)), as.integer(J_max_bs), by = 1) #on commence au moins à 2 parce que 1 c'est pas une bonne solution
  I_hat_bs = sort(I_hat_bs[sapply(I_hat_bs,valid_dim_b_splines, degree)]) #select only the valid dimensions
  
  bs_bool = 0
  J_max_ns = compute_J_max(W, Z, Y, degree, bs_bool)
  I_hat_ns = seq(as.integer(0.1 * (log(J_max_ns)^2)), as.integer(J_max_ns), by = 1) #on commence au moins à 2 parce que 1 c'est pas une bonne solution
  I_hat_ns = sort(I_hat_ns[sapply(I_hat_ns,valid_dim_ns, degree)]) #select only the valid dimensions
  
  # filtrer les valeurs possibles : on enlève ceux où gamma est zéro
  I_hat_bs_new = I_hat_bs[-zero_indices_bs]
  I_hat_ns_new = I_hat_ns[-zero_indices_ns]
  
  if (length(I_hat_bs) > 0){
    c_0 <- 10
    J_opt_lepski_bs <- lepski_chen(I_hat_bs_new, c_0, m_m, list_g_hat_J_bs, W, Z, Y, degree, 
                                   valid_dim_b_splines, create_dyadic_P_splines_bs)
    #j_index_ <- which(J_CV_bs == J_opt_lepski_bs)
    #g_hat_J_bs_lepski = list_g_hat_J_bs[[j_index_]]
    
    # lepski boot
    J_opt_lepskiboot_bs <- lepski_bootstrap(I_hat_bs_new, n_boot_lepski, list_gamma_bs, list_M_boot_bs, 
                                            list_g_hat_J_bs, valid_dim_b_splines, x_evaluation, 
                                            W, Z, Y, degree, create_dyadic_P_splines_bs)
    #j_index_ <- which(J_CV_bs == J_opt_lepskiboot_bs)
    #g_hat_J_bs_lepski_boot = list_g_hat_J_bs[[j_index_]]}
  }
  else{
    J_opt_lepski_bs <- 0
    #g_hat_J_bs_lepski <- rep(0, length(x_evaluation))
    J_opt_lepskiboot_bs <- 0
    #g_hat_J_bs_lepski_boot <- rep(0, length(x_evaluation))
  }
  
  
  if (length(I_hat_ns)>0){
    J_opt_lepski_ns <- lepski_chen(I_hat_ns_new, c_0, m_m, list_g_hat_J_ns, W, Z, Y, degree, 
                                   valid_dim_ns, create_dyadic_P_splines_ns)
    #j_index_ <- which(J_CV_ns == J_opt_lepski_ns)
    #g_hat_J_ns_lepski = list_g_hat_J_ns[[j_index_]]
    
    J_opt_lepskiboot_ns <- lepski_bootstrap(I_hat_ns_new, n_boot_lepski, list_gamma_ns, list_M_boot_ns,
                                            list_g_hat_J_ns, valid_dim_ns, x_evaluation, 
                                            W, Z, Y, degree, create_dyadic_P_splines_ns)

    #j_index_ <- which(J_CV_ns == J_opt_lepskiboot_ns)
    #g_hat_J_ns_lepski_boot = list_g_hat_J_ns[[j_index_]]}
  }
  else{
    J_opt_lepski_ns <- 0
    #g_hat_J_ns_lepski <- rep(0, length(x_evaluation))
    J_opt_lepskiboot_ns <- 0
    #g_hat_J_ns_lepski_boot <- rep(0, length(x_evaluation))
  }
  
  rm(list_M_boot_bs, list_M_boot_ns, W , Y, Z, I_hat_bs, I_hat_ns)
  gc()
  
  
  #return the results
  return(list(J_opt_CV_M_bs = J_opt_CV_M_bs, J_opt_CV_M_ns = J_opt_CV_M_ns,
              J_opt_CVMSE_bs = J_opt_CVMSE_bs, J_opt_CVMSE_ns = J_opt_CVMSE_ns, 
              J_opt_lepski_bs = J_opt_lepski_bs, J_opt_lepski_ns = J_opt_lepski_ns,
              J_opt_lepskiboot_bs = J_opt_lepskiboot_bs, J_opt_lepskiboot_ns = J_opt_lepskiboot_ns,
              list_gamma_bs = list_gamma_bs, list_gamma_ns = list_gamma_ns,
              list_g_hat_J_bs = list_g_hat_J_bs, list_g_hat_J_ns = list_g_hat_J_ns))
              #g_hat_J_bs_CVM = g_hat_J_bs_CVM, g_hat_J_ns_CVM = g_hat_J_ns_CVM,
              #g_hat_J_bs_CVMSE = g_hat_J_bs_CVMSE, g_hat_J_ns_CVMSE = g_hat_J_ns_CVMSE,
              #g_hat_J_bs_lepski = g_hat_J_bs_lepski, g_hat_J_ns_lepski = g_hat_J_ns_lepski,
              #g_hat_J_bs_lepski_boot = g_hat_J_bs_lepski_boot, g_hat_J_ns_lepski_boot = g_hat_J_ns_lepski_boot))
  
}


#### MC ###
MC_j_opt <- function(n_MC, data_param, case, J_CV_bs, J_CV_ns, p_train, x_evaluation,
                     degree, n_boot_lepski){
  
  n_eval = length(x_evaluation)
  list_J_opt_CV_M_bs <- numeric(n_MC)
  list_J_opt_CV_M_ns <- numeric(n_MC)
  list_J_opt_CVMSE_bs <- numeric(n_MC)
  list_J_opt_CVMSE_ns <- numeric(n_MC)
  list_J_opt_lepski_bs <- numeric(n_MC)
  list_J_opt_lepski_ns <- numeric(n_MC)
  list_J_opt_lepskiboot_bs <- numeric(n_MC)
  list_J_opt_lepskiboot_ns <- numeric(n_MC)
  
  #list_J_opt <- vector("list", 8)
  
  list_all_gamma_bs <- vector("list", n_MC)
  list_all_gamma_ns <- vector("list", n_MC)
  
  #list_g_hat_J_bs_CVM <- matrix(NA, nrow = n_MC, ncol = n_eval)
  #list_g_hat_J_ns_CVM <- matrix(NA, nrow = n_MC, ncol = n_eval)
  #list_g_hat_J_bs_CVMSE <- matrix(NA, nrow = n_MC, ncol = n_eval)
  #list_g_hat_J_ns_CVMSE <- matrix(NA, nrow = n_MC, ncol = n_eval)
  #list_g_hat_J_bs_lepski <- matrix(NA, nrow = n_MC, ncol = n_eval)
  #list_g_hat_J_ns_lepski <- matrix(NA, nrow = n_MC, ncol = n_eval)
  #list_g_hat_J_bs_lepski_boot <- matrix(NA, nrow = n_MC, ncol = n_eval)
  #list_g_hat_J_ns_lepski_boot <- matrix(NA, nrow = n_MC, ncol = n_eval)
  
  list_g_hat_J_bs_all <- vector("list", n_MC)
  list_g_hat_J_ns_all <- vector("list", n_MC)
  
  for (n in 1:n_MC){
    print(n)
    simul <- simulate_data_3(data_param, g_sim_3, case)
    res <- J_opt_data_fixed(simul, J_CV_bs, J_CV_ns, p_train, x_evaluation,
                            degree, n_boot_lepski)
    
    #res <- list(J_opt_CV_M_bs = sample(J_CV_bs, 1), J_opt_CV_M_ns = sample(J_CV_bs, 1),
    #            J_opt_CVMSE_bs = sample(J_CV_bs, 1), J_opt_CVMSE_ns = sample(J_CV_bs,1),
    #            J_opt_lepski_bs = sample(J_CV_bs, 1),J_opt_lepski_ns = sample(J_CV_bs, 1),
    #            J_opt_lepskiboot_bs = sample(J_CV_bs, 1), J_opt_lepskiboot_ns = sample(J_CV_bs, 1),
    #            g_hat_J_bs_CVM = rnorm(100, 0, 1), g_hat_J_ns_CVM = rnorm(100, 0, 1), 
    #            g_hat_J_bs_CVMSE = rnorm(100, 0, 1), g_hat_J_ns_CVMSE = rnorm(100, 0, 1),
    #            g_hat_J_bs_lepski = rnorm(100, 0, 1), g_hat_J_ns_lepski =rnorm(100, 0, 1) , 
    #            g_hat_J_bs_lepski_boot= rnorm(100, 0, 1), g_hat_J_ns_lepski_boot= rnorm(100, 0, 1))
                
    list_J_opt_CV_M_bs[n] <- res$J_opt_CV_M_bs
    list_J_opt_CV_M_ns[n] <- res$J_opt_CV_M_ns
    list_J_opt_CVMSE_bs[n] <- res$J_opt_CVMSE_bs
    list_J_opt_CVMSE_ns[n] <- res$J_opt_CVMSE_ns
    list_J_opt_lepski_bs[n] <- res$J_opt_lepski_bs
    list_J_opt_lepski_ns[n] <- res$J_opt_lepski_ns
    list_J_opt_lepskiboot_bs[n] <- res$J_opt_lepskiboot_bs
    list_J_opt_lepskiboot_ns[n] <- res$J_opt_lepskiboot_ns
    
    list_all_gamma_bs[[n]] <- res$list_gamma_bs
    list_all_gamma_ns[[n]] <- res$list_gamma_ns
    
    #list_g_hat_J_bs_CVM[n,] <- res$g_hat_J_bs_CVM
    #list_g_hat_J_ns_CVM[n,] <- res$g_hat_J_ns_CVM
    #list_g_hat_J_bs_CVMSE[n,] <- res$g_hat_J_bs_CVMSE
    #list_g_hat_J_ns_CVMSE[n,] <- res$g_hat_J_ns_CVMSE
    #list_g_hat_J_bs_lepski[n,] <- res$g_hat_J_bs_lepski
    #list_g_hat_J_ns_lepski[n,] <- res$g_hat_J_ns_lepski
    #list_g_hat_J_bs_lepski_boot[n,] <- res$g_hat_J_bs_lepski_boot
    #list_g_hat_J_ns_lepski_boot[n,] <- res$g_hat_J_ns_lepski_boot
    
    list_g_hat_J_bs_all[[n]] <- res$list_g_hat_J_bs
    list_g_hat_J_ns_all[[n]] <- res$list_g_hat_J_ns
    
    #rm(simul, res)
    #gc()
  }
  
  return(list(list_J_opt_CV_M_bs = list_J_opt_CV_M_bs, list_J_opt_CV_M_ns = list_J_opt_CV_M_ns,
         list_J_opt_CVMSE_bs = list_J_opt_CVMSE_bs, list_J_opt_CVMSE_ns = list_J_opt_CVMSE_ns,
         list_J_opt_lepski_bs = list_J_opt_lepski_bs, list_J_opt_lepski_ns = list_J_opt_lepski_ns,
         list_J_opt_lepskiboot_bs = list_J_opt_lepskiboot_bs, list_J_opt_lepskiboot_ns = list_J_opt_lepskiboot_ns,
         list_g_hat_J_bs_all = list_g_hat_J_bs_all,list_g_hat_J_ns_all = list_g_hat_J_ns_all
         #list_all_gamma_bs = list_all_gamma_bs, list_all_gamma_ns = list_all_gamma_ns,
         #list_g_hat_J_bs_CVM = list_g_hat_J_bs_CVM, list_g_hat_J_ns_CVM = list_g_hat_J_ns_CVM,
         #list_g_hat_J_bs_CVMSE = list_g_hat_J_bs_CVMSE, list_g_hat_J_ns_CVMSE = list_g_hat_J_ns_CVMSE,
         #list_g_hat_J_bs_lepski = list_g_hat_J_bs_lepski, list_g_hat_J_ns_lepski = list_g_hat_J_ns_lepski,
         #list_g_hat_J_bs_lepski_boot = list_g_hat_J_bs_lepski_boot, list_g_hat_J_ns_lepski_boot = list_g_hat_J_ns_lepski_boot))
  ))
         }



J_bs_ <- c(5, 7, 11, 19, 35)
J_ns_ <- c(3, 5, 9, 17, 33)
p = 0.5
x_eval = seq(-2, 2, length.out = 100)
deg = 3
n_boot = 10
simul <- simulate_data_3(c(200, 0.5, 0.9), g_sim_3, 2)

test <- J_opt_data_fixed(simul, J_bs_, J_ns_, p, x_eval, deg, n_boot)


test <- MC_j_opt(1, c(200, 0.5, 0.9), 2, J_bs_, J_ns_, p, x_eval,
         deg, n_boot)


