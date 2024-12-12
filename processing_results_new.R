#### Packages ####
#setwd("C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/simulation_results/3")
setwd("C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/github/master_thesis")
source("source_file_common_all.R")
library(splines)
library(MASS)
library(caret)
library(expm)
library(foreach)
library(ggplot2)
library(dplyr)
library(stringr)
library(tibble)
library(xtable)
library(kableExtra)

setwd("C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/simulation_results/3")

#### Load results ####
#load("opt_5_degree3_ptrain0.5_nboot100_rhozw0.9_rhouv0.5_case2_n1000.R")
file_list <- list.files(path = "C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/simulation_results/3",
                        pattern = "^opt_.*\\.R$", full.names = FALSE)

all_data <- lapply(file_list, function(filepath) {
  # Extract the filename without the full path
  filename <- basename(filepath)
  name <- str_remove(filename, "\\.R$")  # Remove the ".R" extension
  
  # Load the data from the file
  load(filepath)  # Assumes "results" object is loaded
  
  # Return the results, named by the filename
  list(name = name, data = results)
})

all_data <- setNames(lapply(all_data, `[[`, "data"), sapply(all_data, `[[`, "name"))

#### Functions to process the data ####
#compute missing iterations number (except for each method : done later)
missing_number <- function(n_MC, degree, p_train, n_boot, rhozw, rhouv, case, n_values){
  simu_name = paste("opt_", n_MC , "_degree", degree, "_ptrain", p_train, "_nboot",
                    n_boot, "_rhozw" , rhozw,"_rhouv", rhouv , "_case", case,
                    "_n", n_values ,sep = "") 
  
  is_null <- sapply(all_data[[simu_name]], is.null)
  n_null <- sum(is_null) #first the buggs so it did not even run : so did not run for any method
  all_data[[simu_name]] <- all_data[[simu_name]][!is_null] #remove them
  n_MC = length(all_data[[simu_name]])
  
  #gammas
  n_gamma_bs_5 <- n_null
  n_gamma_bs_7 <- n_null
  n_gamma_bs_11 <- n_null
  n_gamma_bs_19 <- n_null
  n_gamma_bs_35 <- n_null
  
  n_gamma_ns_3 <- n_null
  n_gamma_ns_5 <- n_null
  n_gamma_ns_9 <- n_null
  n_gamma_ns_17 <- n_null
  n_gamma_ns_33 <- n_null
  
  for (n in 1:n_MC){
    list_gamma_bs <- all_data[[simu_name]][[n]]$list_gamma_bs
    if (sum(list_gamma_bs[[1]])==0){
      n_gamma_bs_5 = n_gamma_bs_5 + 1
    }
    if (sum(list_gamma_bs[[2]])==0){
      n_gamma_bs_7 = n_gamma_bs_7 + 1
    }
    if (sum(list_gamma_bs[[3]])==0){
      n_gamma_bs_11 = n_gamma_bs_11 + 1
    }
    if (sum(list_gamma_bs[[4]])==0){
      n_gamma_bs_19 = n_gamma_bs_19 + 1
    }
    if (sum(list_gamma_bs[[5]])==0){
      n_gamma_bs_35 = n_gamma_bs_35 +1 
    }
    
    list_gamma_ns <- all_data[[simu_name]][[n]]$list_gamma_ns
    if (sum(list_gamma_ns[[1]])==0){
      n_gamma_ns_3 = n_gamma_ns_3 + 1
    }
    if (sum(list_gamma_ns[[2]])==0){
      n_gamma_ns_5 = n_gamma_ns_5 + 1
    }
    if (sum(list_gamma_ns[[3]])==0){
      n_gamma_ns_9 = n_gamma_ns_9 + 1
    }
    if (sum(list_gamma_ns[[4]])==0){
      n_gamma_ns_17 = n_gamma_ns_17 + 1
    }
    if (sum(list_gamma_ns[[5]])==0){
      n_gamma_ns_33 = n_gamma_ns_33 +1 
    }    
    
  }

  
  return(list(n_null = n_null, n_gamma_bs_5 = n_gamma_bs_5, 
              n_gamma_bs_7 = n_gamma_bs_7 , n_gamma_bs_11 = n_gamma_bs_11,
              n_gamma_bs_19 = n_gamma_bs_19, n_gamma_bs_35 = n_gamma_bs_35,
              n_gamma_ns_3 = n_gamma_ns_3, n_gamma_ns_5 = n_gamma_ns_5, 
              n_gamma_ns_9 = n_gamma_ns_9, n_gamma_ns_17 = n_gamma_ns_17,
              n_gamma_ns_33 = n_gamma_ns_33))
}

reshape_simul_all <- function(simul_all){
  n_MC = length(simul_all)
  n_values = length(simul_all[[1]]$Y)
  simul_all_sorted <- list(length=3)
  for (i in 1:2){
    simul_all_sorted$Y <- matrix(0, n_MC, n_values)
    simul_all_sorted$W <- matrix(0, n_MC, n_values)
  }
  for (n in 1:n_MC){
    simul_all_sorted$Y[n, ] <- simul_all[[n]]$Y
    simul_all_sorted$W[n, ] <- simul_all[[n]]$W 
  }
  return(simul_all_sorted)
}



create_list_to_analyze <- function(n_MC, degree, p_train, n_boot, rhozw, rhouv, case, n_values, simul_all){
  simu_name = paste("opt_", n_MC , "_degree", degree, "_ptrain", p_train, "_nboot",
                    n_boot, "_rhozw" , rhozw,"_rhouv", rhouv , "_case", case,
                    "_n", n_values ,sep = "") 
  is_null <- sapply(all_data[[simu_name]], is.null) #length 2000 attention, faudra garder les indices et checker avec la data
  all_data[[simu_name]] <- all_data[[simu_name]][!is_null]
  
  #simul_all
  simul_all <- simul_all[!is_null]
  simul_all_reshaped <- reshape_simul_all(simul_all)
  
  n_MC = length(all_data[[simu_name]])
  
  all_lists <- list()
  all_lists$list_J_opt_CV_M_bs <- numeric(length = n_MC)
  all_lists$list_J_opt_CV_M_ns <- numeric(length = n_MC)
  all_lists$list_J_opt_CVMSE_bs <- numeric(length = n_MC)
  all_lists$list_J_opt_CVMSE_ns <- numeric(length = n_MC)
  all_lists$list_J_opt_lepski_bs <- numeric(length = n_MC)
  all_lists$list_J_opt_lepski_ns <- numeric(length = n_MC)
  all_lists$list_J_opt_lepskiboot_bs <- numeric(length = n_MC)
  all_lists$list_J_opt_lepskiboot_ns <- numeric(length = n_MC)
  all_lists$list_gamma_bs <- vector('list', length = n_MC)
  all_lists$list_gamma_ns <-  vector('list', length = n_MC)
  all_lists$list_g_hat_J_bs <- vector('list', length = n_MC)
  all_lists$list_g_hat_J_ns <-  vector('list', length = n_MC)
  all_lists$list_matrix_g_hat_J_bs <- vector("list", length = 5)
  all_lists$list_matrix_g_hat_J_ns <- vector("list", length = 5)
  all_lists$matrix_W <- simul_all_reshaped$W
  all_lists$matrix_Y <-  simul_all_reshaped$Y
  #all_lists$matrix_Z <-  matrix(0, nrow = n_MC, ncol = n_values)
  
  for (n in 1:n_MC){
    all_lists$list_J_opt_CV_M_bs[n] <- all_data[[simu_name]][[n]]$J_opt_CV_M_bs
    all_lists$list_J_opt_CV_M_ns[n] <- all_data[[simu_name]][[n]]$J_opt_CV_M_ns
    all_lists$list_J_opt_CVMSE_bs[n] <- all_data[[simu_name]][[n]]$J_opt_CVMSE_bs
    all_lists$list_J_opt_CVMSE_ns[n] <- all_data[[simu_name]][[n]]$J_opt_CVMSE_ns
    all_lists$list_J_opt_lepski_bs[n] <- all_data[[simu_name]][[n]]$J_opt_lepski_bs
    all_lists$list_J_opt_lepski_ns[n] <- all_data[[simu_name]][[n]]$J_opt_lepski_ns
    all_lists$list_J_opt_lepskiboot_bs[n] <- all_data[[simu_name]][[n]]$J_opt_lepskiboot_bs
    all_lists$list_J_opt_lepskiboot_ns[n] <- all_data[[simu_name]][[n]]$J_opt_lepskiboot_ns
    all_lists$list_gamma_bs[[n]] <- all_data[[simu_name]][[n]]$list_gamma_bs
    all_lists$list_gamma_ns[[n]] <- all_data[[simu_name]][[n]]$list_gamma_ns
    all_lists$list_g_hat_J_bs[[n]] <- all_data[[simu_name]][[n]]$list_g_hat_J_bs
    all_lists$list_g_hat_J_ns[[n]] <- all_data[[simu_name]][[n]]$list_g_hat_J_ns
  }
  
  for (i in 1:5){
    all_lists$list_matrix_g_hat_J_bs[[i]] <- matrix(0, nrow = n_MC, ncol = n_eval)
    all_lists$list_matrix_g_hat_J_ns[[i]] <- matrix(0, nrow = n_MC, ncol = n_eval)
    for (n in 1:n_MC){
      all_lists$list_matrix_g_hat_J_bs[[i]][n,] <- all_data[[simu_name]][[n]]$list_g_hat_J_bs[[i]]
      all_lists$list_matrix_g_hat_J_ns[[i]][n, ] <- all_data[[simu_name]][[n]]$list_g_hat_J_ns[[i]]
      }
  }
  
  # create the lists for each method of values on x_eval obtained
  all_lists$list_est_values_CV_M_bs <- matrix(0, nrow = n_MC, ncol = n_eval)
  all_lists$list_est_values_CV_M_ns <- matrix(0, nrow = n_MC, ncol = n_eval)
  all_lists$list_est_values_CVMSE_bs <- matrix(0, nrow = n_MC, ncol = n_eval)
  all_lists$list_est_values_CVMSE_ns <- matrix(0, nrow = n_MC, ncol = n_eval)
  all_lists$list_est_values_lepski_bs <- matrix(0, nrow = n_MC, ncol = n_eval)
  all_lists$list_est_values_lepski_ns <- matrix(0, nrow = n_MC, ncol = n_eval)
  all_lists$list_est_values_lepskiboot_bs <- matrix(0, nrow = n_MC, ncol = n_eval)
  all_lists$list_est_values_lepskiboot_ns <- matrix(0, nrow = n_MC, ncol = n_eval)

  #check for NULL values 
  for (n in 1:n_MC){
    #bs
    index = which(J_bs == all_lists$list_J_opt_CV_M_bs[n])
    if (length(index)>0){
      all_lists$list_est_values_CV_M_bs[n,] <- all_lists$list_g_hat_J_bs[[n]][index][[1]]}
    
    index = which(J_bs == all_lists$list_J_opt_CVMSE_bs[n])
    if (length(index)>0){
      all_lists$list_est_values_CVMSE_bs[n,] <- all_lists$list_g_hat_J_bs[[n]][index][[1]]}
    
    index = which(J_bs == all_lists$list_J_opt_lepski_bs[n])
    if (length(index)>0){
      all_lists$list_est_values_lepski_bs[n,] <- all_lists$list_g_hat_J_bs[[n]][index][[1]]}

    
    index = which(J_bs == all_lists$list_J_opt_lepskiboot_bs[n])
    if (length(index)>0){
      all_lists$list_est_values_lepskiboot_bs[n,] <- all_lists$list_g_hat_J_bs[[n]][index][[1]]}

    
    index = which(J_ns == all_lists$list_J_opt_CV_M_ns[n])
    if (length(index)>0){
      all_lists$list_est_values_CV_M_ns[n,] <- all_lists$list_g_hat_J_ns[[n]][index][[1]]}

    index = which(J_ns == all_lists$list_J_opt_CVMSE_ns[n])
    if (length(index)>0){
      all_lists$list_est_values_CVMSE_ns[n,] <- all_lists$list_g_hat_J_ns[[n]][index][[1]]}
    
    index = which(J_ns == all_lists$list_J_opt_lepski_ns[n])
    if (length(index)>0){
      all_lists$list_est_values_lepski_ns[n,] <- all_lists$list_g_hat_J_ns[[n]][index][[1]]}

    index = which(J_ns == all_lists$list_J_opt_lepskiboot_ns[n])
    if (length(index)>0){
      all_lists$list_est_values_lepskiboot_ns[n,] <- all_lists$list_g_hat_J_ns[[n]][index][[1]]}
    }
  
  # remove unnecessary stuff
  #names_to_remove <- c("list_gamma_bs", "list_gamma_ns", "list_g_hat_J_bs", "list_g_hat_J_ns")
  #all_lists <- all_lists[setdiff(names(all_lists), names_to_remove)]
  
  return(all_lists)
}




# filter the observations that equal 0
filter_res <- function(res_to_analyze){
  #CVM bs
  res_CVM_bs <- list()
  zero_indices <- which(res_to_analyze$list_J_opt_CV_M_bs == 0) #les J qui valent 0 
  n_CVM_bs <- length(zero_indices)
  if (n_CVM_bs>0){
    res_CVM_bs$list_J_opt_CV_M_bs <- res_to_analyze$list_J_opt_CV_M_bs[-zero_indices]
    res_CVM_bs$list_est_values_CV_M_bs <- res_to_analyze$list_est_values_CV_M_bs[-zero_indices,]
    # res_CVM_bs$matrix_W <- res_to_analyze$matrix_W[-zero_indices,]
    # res_CVM_bs$matrix_Z <- res_to_analyze$matrix_Z[-zero_indices,]
    # res_CVM_bs$matrix_Y <- res_to_analyze$matrix_Y[-zero_indices,]
    }
  else{
    res_CVM_bs$list_J_opt_CV_M_bs <- res_to_analyze$list_J_opt_CV_M_bs
    res_CVM_bs$list_est_values_CV_M_bs <- res_to_analyze$list_est_values_CV_M_bs
    # res_CVM_bs$matrix_W <- res_to_analyze$matrix_W
    # res_CVM_bs$matrix_Z <- res_to_analyze$matrix_Z
    # res_CVM_bs$matrix_Y <- res_to_analyze$matrix_Y
  }
  
  #CVM ns
  res_CVM_ns <- list()
  zero_indices <- which(res_to_analyze$list_J_opt_CV_M_ns == 0) #les J qui valent 0 
  n_CVM_ns <- length(zero_indices)
  if (n_CVM_ns>0){
    res_CVM_ns$list_J_opt_CV_M_ns <- res_to_analyze$list_J_opt_CV_M_ns[-zero_indices]
    res_CVM_ns$list_est_values_CV_M_ns <- res_to_analyze$list_est_values_CV_M_ns[-zero_indices,]
    res_CVM_ns$matrix_W <- res_to_analyze$matrix_W[-zero_indices,]
    res_CVM_ns$matrix_Z <- res_to_analyze$matrix_Z[-zero_indices,]
    res_CVM_ns$matrix_Y <- res_to_analyze$matrix_Y[-zero_indices,]
    }
  else{
    res_CVM_ns$list_J_opt_CV_M_ns <- res_to_analyze$list_J_opt_CV_M_ns
    res_CVM_ns$list_est_values_CV_M_ns <- res_to_analyze$list_est_values_CV_M_ns
    res_CVM_ns$matrix_W <- res_to_analyze$matrix_W
    res_CVM_ns$matrix_Z <- res_to_analyze$matrix_Z
    res_CVM_ns$matrix_Y <- res_to_analyze$matrix_Y
  }
  
  #CVMSE bs
  res_CVMSE_bs <- list()
  zero_indices <- which(res_to_analyze$list_J_opt_CVMSE_bs == 0) #les J qui valent 0 
  n_CVMSE_bs <- length(zero_indices)
  if (n_CVMSE_bs>0){
    res_CVMSE_bs$list_J_opt_CVMSE_bs <- res_to_analyze$list_J_opt_CVMSE_bs[-zero_indices]
    res_CVMSE_bs$list_est_values_CVMSE_bs <- res_to_analyze$list_est_values_CVMSE_bs[-zero_indices,]
    res_CVMSE_bs$matrix_W <- res_to_analyze$matrix_W[-zero_indices,]
    res_CVMSE_bs$matrix_Z <- res_to_analyze$matrix_Z[-zero_indices,]
    res_CVMSE_bs$matrix_Y <- res_to_analyze$matrix_Y[-zero_indices,]
    }
  else{
    res_CVMSE_bs$list_J_opt_CVMSE_bs <- res_to_analyze$list_J_opt_CVMSE_bs
    res_CVMSE_bs$list_est_values_CVMSE_bs <- res_to_analyze$list_est_values_CVMSE_bs
    res_CVMSE_bs$matrix_W <- res_to_analyze$matrix_W
    res_CVMSE_bs$matrix_Z <- res_to_analyze$matrix_Z
    res_CVMSE_bs$matrix_Y <- res_to_analyze$matrix_Y
  }
  
  #CVMSE ns
  res_CVMSE_ns <- list()
  zero_indices <- which(res_to_analyze$list_J_opt_CVMSE_ns == 0) #les J qui valent 0 
  n_CVMSE_ns <- length(zero_indices)
  if (n_CVMSE_ns>0){
    res_CVMSE_ns$list_J_opt_CVMSE_ns <- res_to_analyze$list_J_opt_CVMSE_ns[-zero_indices]
    res_CVMSE_ns$list_est_values_CVMSE_ns <- res_to_analyze$list_est_values_CVMSE_ns[-zero_indices,]
    res_CVMSE_ns$matrix_W <- res_to_analyze$matrix_W[-zero_indices,]
    res_CVMSE_ns$matrix_Z <- res_to_analyze$matrix_Z[-zero_indices,]
    res_CVMSE_ns$matrix_Y <- res_to_analyze$matrix_Y[-zero_indices,]
    }
  else{
    res_CVMSE_ns$list_J_opt_CVMSE_ns <- res_to_analyze$list_J_opt_CVMSE_ns
    res_CVMSE_ns$list_est_values_CVMSE_ns <- res_to_analyze$list_est_values_CVMSE_ns
    res_CVMSE_ns$matrix_W <- res_to_analyze$matrix_W
    res_CVMSE_ns$matrix_Z <- res_to_analyze$matrix_Z
    res_CVMSE_ns$matrix_Y <- res_to_analyze$matrix_Y
    }
  
  #lepski bs
  res_lepski_bs <- list()
  zero_indices <- which(res_to_analyze$list_J_opt_lepski_bs == 0) #les J qui valent 0 
  n_lepski_bs <- length(zero_indices)
  if (n_lepski_bs>0){
    res_lepski_bs$list_J_opt_lepski_bs <- res_to_analyze$list_J_opt_lepski_bs[-zero_indices]
    res_lepski_bs$list_est_values_lepski_bs <- res_to_analyze$list_est_values_lepski_bs[-zero_indices,]
    res_lepski_bs$matrix_W <- res_to_analyze$matrix_W[-zero_indices,]
    res_lepski_bs$matrix_Z <- res_to_analyze$matrix_Z[-zero_indices,]
    res_lepski_bs$matrix_Y <- res_to_analyze$matrix_Y[-zero_indices,]
    }
  else{
    res_lepski_bs$list_J_opt_lepski_bs <- res_to_analyze$list_J_opt_lepski_bs
    res_lepski_bs$list_est_values_lepski_bs <- res_to_analyze$list_est_values_lepski_bs
    res_lepski_bs$matrix_W <- res_to_analyze$matrix_W
    res_lepski_bs$matrix_Z <- res_to_analyze$matrix_Z
    res_lepski_bs$matrix_Y <- res_to_analyze$matrix_Y
    }

  #lepski ns
  res_lepski_ns <- list()
  zero_indices <- which(res_to_analyze$list_J_opt_lepski_ns == 0) #les J qui valent 0 
  n_lepski_ns <- length(zero_indices)
  if (n_lepski_ns>0){
    res_lepski_ns$list_J_opt_lepski_ns <- res_to_analyze$list_J_opt_lepski_ns[-zero_indices]
    res_lepski_ns$list_est_values_lepski_ns <- res_to_analyze$list_est_values_lepski_ns[-zero_indices,]
    res_lepski_ns$matrix_W <- res_to_analyze$matrix_W[-zero_indices,]
    res_lepski_ns$matrix_Z <- res_to_analyze$matrix_Z[-zero_indices,]
    res_lepski_ns$matrix_Y <- res_to_analyze$matrix_Y[-zero_indices,]
    }
  else{
    res_lepski_ns$list_J_opt_lepski_ns <- res_to_analyze$list_J_opt_lepski_ns
    res_lepski_ns$list_est_values_lepski_ns <- res_to_analyze$list_est_values_lepski_ns
    res_lepski_ns$matrix_W <- res_to_analyze$matrix_W
    res_lepski_ns$matrix_Z <- res_to_analyze$matrix_Z
    res_lepski_ns$matrix_Y <- res_to_analyze$matrix_Y
    }  
  
  #lepskiboot bs
  res_lepskiboot_bs <- list()
  zero_indices <- which(res_to_analyze$list_J_opt_lepskiboot_bs == 0) #les J qui valent 0 
  n_lepskiboot_bs <- length(zero_indices)
  if (n_lepskiboot_bs>0){
    res_lepskiboot_bs$list_J_opt_lepskiboot_bs <- res_to_analyze$list_J_opt_lepskiboot_bs[-zero_indices]
    res_lepskiboot_bs$list_est_values_lepskiboot_bs <- res_to_analyze$list_est_values_lepskiboot_bs[-zero_indices,]
    res_lepskiboot_bs$matrix_W <- res_to_analyze$matrix_W[-zero_indices,]
    res_lepskiboot_bs$matrix_Z <- res_to_analyze$matrix_Z[-zero_indices,]
    res_lepskiboot_bs$matrix_Y <- res_to_analyze$matrix_Y[-zero_indices,]
    }
  else{
    res_lepskiboot_bs$list_J_opt_lepskiboot_bs <- res_to_analyze$list_J_opt_lepskiboot_bs
    res_lepskiboot_bs$list_est_values_lepskiboot_bs <- res_to_analyze$list_est_values_lepskiboot_bs
    res_lepskiboot_bs$matrix_W <- res_to_analyze$matrix_W
    res_lepskiboot_bs$matrix_Z <- res_to_analyze$matrix_Z
    res_lepskiboot_bs$matrix_Y <- res_to_analyze$matrix_Y
    }
  
  #lepskiboot ns
  res_lepskiboot_ns <- list()
  zero_indices <- which(res_to_analyze$list_J_opt_lepskiboot_ns == 0) #les J qui valent 0 
  n_lepskiboot_ns <- length(zero_indices)
  if (n_lepskiboot_ns>0){
    res_lepskiboot_ns$list_J_opt_lepskiboot_ns <- res_to_analyze$list_J_opt_lepskiboot_ns[-zero_indices]
    res_lepskiboot_ns$list_est_values_lepskiboot_ns <- res_to_analyze$list_est_values_lepskiboot_ns[-zero_indices,]
    res_lepskiboot_ns$matrix_W <- res_to_analyze$matrix_W[-zero_indices,]
    res_lepskiboot_ns$matrix_Z <- res_to_analyze$matrix_Z[-zero_indices,]
    res_lepskiboot_ns$matrix_Y <- res_to_analyze$matrix_Y[-zero_indices,]
    }
  else{
    res_lepskiboot_ns$list_J_opt_lepskiboot_ns <- res_to_analyze$list_J_opt_lepskiboot_ns
    res_lepskiboot_ns$list_est_values_lepskiboot_ns <- res_to_analyze$list_est_values_lepskiboot_ns
    res_lepskiboot_ns$matrix_W <- res_to_analyze$matrix_W
    res_lepskiboot_ns$matrix_Z <- res_to_analyze$matrix_Z
    res_lepskiboot_ns$matrix_Y <- res_to_analyze$matrix_Y
    }
  
  return(list(res_CVM_bs = res_CVM_bs, res_CVM_ns = res_CVM_ns,
              res_CVMSE_bs = res_CVMSE_bs, res_CVMSE_ns = res_CVMSE_ns,
              res_lepski_bs = res_lepski_bs, res_lepski_ns = res_lepski_ns,
              res_lepskiboot_bs = res_lepskiboot_bs, res_lepskiboot_ns = res_lepskiboot_ns,
              n_CVM_bs = n_CVM_bs, n_CVM_ns = n_CVM_ns, n_CVMSE_bs = n_CVMSE_bs,
              n_CVMSE_ns = n_CVMSE_ns, n_lepski_bs = n_lepski_bs, n_lepski_ns = n_lepski_ns,
              n_lepskiboot_bs = n_lepskiboot_bs, n_lepskiboot_ns = n_lepskiboot_ns
              ))
  
}





filter_matrix_pred <- function(matrix_pred) {
  # Check which rows have a sum not equal to 0
  non_zero_rows <- apply(matrix_pred, 1, function(row) sum(row) != 0)
  
  # Subset the matrix to include only the non-zero rows
  new_matrix <- matrix_pred[non_zero_rows, , drop = FALSE]
  
  return(new_matrix)
}


# Compute performances
compute_all_perf <- function(mat_g_hat_on_x, g_0_on_x, data_algo){
  n_eval = length(g_0_on_x)
  n_MC = nrow(mat_g_hat_on_x)

  #MSE
  MSE = 0 
  for (n in 1:n_MC){
    MSE = MSE + sum((mat_g_hat_on_x[n,] - g_0_on_x)^{2})/n_eval
  }
  MSE = MSE/n_MC
  
  #var
  var = 0
  n_eval = length(mat_g_hat_on_x[1,])
  avg <- rep(0, n_eval)
  for (x in 1:n_eval){
    for (n in 1:n_MC){
      avg[x] <- avg[x] + mat_g_hat_on_x[n,][x]/n_MC
    }}
  for (n in 1:n_MC){
    var = var + sum((mat_g_hat_on_x[n,] - avg)^{2})/n_eval
  }
  var = var/(n_MC)
  
  
  #bias 
  bias = 0 
  for (n in 1:n_MC){
    bias = bias + sum((g_0_on_x - avg)^{2})/n_eval
  }
  bias = bias/(n_MC)
  
  #supnorm
  sup_norm_vect <- rep(0, n_MC)
  for (n in 1:n_MC){
    sup_norm_vect[n] = max(abs(mat_g_hat_on_x[n,] - g_0_on_x))
  }
  sup_norm = mean(sup_norm_vect)

  #M A TESTER
  M_vect <- rep(0, n_MC)
  matrix_Y <- data_algo$matrix_Y
  matrix_Z <- data_algo$matrix_Z
  matrix_W <- data_algo$matrix_W
  for (n in 1:n_MC){
    Omega <- create_W(matrix_W[n,])
    M_vect[n] <- calcul_M_g_hat_test_sample(mat_g_hat_on_x[n,], Omega, n_eval, matrix_Y[n,])
  }
  M = mean(M_vect)
  #list(MSE = MSE, var = var, bias = bias, supnorm = sup_norm, M = M)
  return(c(M, sup_norm, MSE, bias, var))
}


compute_perf <- function(res_filtered, g_0_on_x){
  measures_CVM_bs <- compute_all_perf(res_filtered[["res_CVM_bs"]]$list_est_values_CV_M_bs, g_0_on_x, res_filtered[["res_CVM_bs"]])
  measures_CVM_ns <- compute_all_perf(res_filtered[["res_CVM_ns"]]$list_est_values_CV_M_ns, g_0_on_x, res_filtered[["res_CVM_ns"]])
  measures_CVMSE_bs <- compute_all_perf(res_filtered[["res_CVMSE_bs"]]$list_est_values_CVMSE_bs, g_0_on_x, res_filtered[["res_CVMSE_bs"]])
  measures_CVMSE_ns <- compute_all_perf(res_filtered[["res_CVMSE_ns"]]$list_est_values_CVMSE_ns, g_0_on_x, res_filtered[["res_CVMSE_ns"]])
  measures_lepski_bs <- compute_all_perf(res_filtered[["res_lepski_bs"]]$list_est_values_lepski_bs, g_0_on_x, res_filtered[["res_lepski_bs"]])
  measures_lepski_ns <- compute_all_perf(res_filtered[["res_lepski_ns"]]$list_est_values_lepski_ns, g_0_on_x, res_filtered[["res_lepski_ns"]])
  measures_lepskiboot_bs <- compute_all_perf(res_filtered[["res_lepskiboot_bs"]]$list_est_values_lepskiboot_bs, g_0_on_x, res_filtered[["res_lepskiboot_bs"]])
  measures_lepskiboot_ns <- compute_all_perf(res_filtered[["res_lepskiboot_ns"]]$list_est_values_lepskiboot_ns, g_0_on_x, res_filtered[["res_lepskiboot_ns"]])
  
  return(list(measures_CVM_bs = measures_CVM_bs, measures_CVM_ns = measures_CVM_ns,
         measures_CVMSE_bs = measures_CVMSE_bs, measures_CVMSE_ns = measures_CVMSE_ns,
         measures_lepski_bs = measures_lepski_bs, measures_lepski_ns = measures_lepski_ns,
         measures_lepskiboot_bs = measures_lepskiboot_bs, measures_lepskiboot_ns= measures_lepskiboot_ns))
  
}


#perf J
compute_perf_J_sub <- function(g_eval, g_on_x, avg){#rajouter le Y et le W
  MSE = mean((g_eval - g_on_x)^{2})
  var = mean((g_eval - avg)^{2})
  bias = mean((g_on_x - avg)^{2})
  sup_norm = max(abs(g_eval - g_on_x))
  M = 0 # TO DO à réfléchir comment récupérer les données : 
  
  return(c(M, sup_norm, MSE, bias, var))
}

compute_perf_J <- function(value_bs_all_J, value_ns_all_J, g_on_x){
  to_return <- list()

  for (i in 1:5){
    matrix_eval <- value_bs_all_J[[i]]
    n_row = nrow(matrix_eval)
    if (n_row > 0){
      avg <- colMeans(matrix_eval)
      matrix_perf <- matrix(0, nrow = n_row, ncol = 5)
      for (n in 1:n_row){
        matrix_perf[n,] <- compute_perf_J_sub(matrix_eval[n,], g_on_x, avg)
      }
      perf <- colMeans(matrix_perf)
      }
    else{
      perf <- rep(999, 5)
    }
    if (i == 1){
      to_return$avg_perf_5_bs <- perf
    }
    if (i == 2){
      to_return$avg_perf_7_bs <- perf
    }
    if (i == 3){
      to_return$avg_perf_11_bs <- perf
    }
    if (i == 4){
      to_return$avg_perf_19_bs <- perf
    }
    if (i == 5){
      to_return$avg_perf_35_bs <- perf
    }}

  for (i in 1:5){
    matrix_eval <- value_ns_all_J[[i]]
    n_row = nrow(matrix_eval)
    if (n_row > 0){
      avg <- colMeans(matrix_eval)
      matrix_perf <- matrix(0, nrow = n_row, ncol = 5)
      for (n in 1:n_row){
        matrix_perf[n,] <- compute_perf_J_sub(matrix_eval[n,], g_on_x, avg)
      }
      perf <- colMeans(matrix_perf)
    }
    else{
      perf <- rep(999, 5)
    }
    if (i == 1){
      to_return$avg_perf_3_ns <- perf
    }
    if (i == 2){
      to_return$avg_perf_5_ns <- perf
    }
    if (i == 3){
      to_return$avg_perf_9_ns <- perf
    }
    if (i == 4){
      to_return$avg_perf_17_ns <- perf
    }
    if (i == 5){
      to_return$avg_perf_33_ns <- perf
    }}
  
  return(to_return)
  
  # return(list(avg_perf_5_bs = avg_perf_5_bs, avg_perf_7_bs = avg_perf_7_bs,
  #             avg_perf_11_bs = avg_perf_11_bs, avg_perf_19_bs = avg_perf_19_bs,
  #             avg_perf_35_bs = avg_perf_35_bs,
  #             avg_perf_3_ns = avg_perf_3_ns, avg_perf_5_ns = avg_perf_5_ns,
  #             avg_perf_9_ns = avg_perf_9_ns, avg_perf_17_ns = avg_perf_17_ns, 
  #             avg_perf_33_ns = avg_perf_33_ns))

}




#### Data process tout d'un coup ####
create_df_measures <- function(n_MC, degree, p_train, n_boot, rhozw, rhouv, case, n_values, g_on_x){
  #load simulation data 
  data_name = paste("data_", n_MC, "_rhozw" , rhozw,"_rhouv", rhouv , "_case", case, "_n", n_values, ".R" ,sep = "")
  simul_all <- get(load(data_name))
  
  #missing data 
  missing_iter <- missing_number(n_MC, degree, p_train, n_boot, rhozw, rhouv, case, n_values)
  n_null = missing_iter$n_null #not even computed
  
  #algorithms, filter 
  res_1 <- create_list_to_analyze(n_MC, degree, p_train, n_boot, rhozw, rhouv, case, n_values, simul_all)
  res_1_filter <- filter_res(res_1)
  
  values_bs_allJ = res_1$list_matrix_g_hat_J_bs
  values_ns_allJ = res_1$list_matrix_g_hat_J_ns
  
  for (i in 1:5){
    values_bs_allJ[[i]] <- filter_matrix_pred(values_bs_allJ[[i]])
    values_ns_allJ[[i]] <- filter_matrix_pred(values_ns_allJ[[i]])
  }
  
  #perf algos
  perf_algos <- compute_perf(res_1_filter, g_on_x)
  
  rm(res_1)
  gc()
  
  #perf J fixed
  #perf_J <- compute_perf_J(res_1, g_on_x) #: plus qu'à faire lui et j'aurai tout 
  perf_J <- compute_perf_J(values_bs_allJ, values_ns_allJ, g_on_x)
  
  # missing iterations on the algos
  n_CVM_bs = res_1_filter$n_CVM_bs + n_null
  n_CVM_ns = res_1_filter$n_CVM_ns + n_null
  n_CVMSE_bs = res_1_filter$n_CVMSE_bs + n_null
  n_CVMSE_ns = res_1_filter$n_CVMSE_ns + n_null
  n_lepski_bs = res_1_filter$n_lepski_bs + n_null
  n_lepski_ns = res_1_filter$n_lepski_ns + n_null
  n_lepskiboot_bs = res_1_filter$n_lepskiboot_bs + n_null
  n_lepskiboot_ns = res_1_filter$n_lepskiboot_ns + n_null
  
  
  vect_missing_algos <- c(n_CVM_bs, n_CVM_ns, n_CVMSE_bs, n_CVMSE_ns,
                          n_lepski_bs, n_lepski_ns, n_lepskiboot_bs, n_lepskiboot_ns)
  
  #une fois que j'ai tout : faire sortir seulement un df avec toutes les performances
  
  columns <- c("M", "sup_norm", "MSE", "bias", "var", "missing")
  rows <- c("J_bs = 5", "J_bs = 7", "J_bs = 11", "J_bs = 19", "J_bs = 35", 
            "J_ns = 3", "J_ns = 5", "J_ns = 9", "J_ns = 17", "J_ns = 33", 
            "CVM_bs", "CVM_ns", "CVMSE_bs", "CVMSE_ns", 
            "lepski_bs", "lepski_ns", "lepskiboot_bs", "lepskiboot_ns")
  df <- data.frame(matrix(ncol = length(columns), nrow = length(rows)))
  colnames(df) <- columns
  rownames(df) <- rows
  
  
  # vector_bs
  values <- 2^(1:5) + 3
  vect_bs <- sapply(values, function(value) {
    eval(parse(text = paste0("missing_iter$n_gamma_bs_", value)))
  })
  values <- 2^(1:5) + 1
  vect_ns <- sapply(values, function(value) {
    eval(parse(text = paste0("missing_iter$n_gamma_ns_", value)))
  })
  missing_vector <- c(vect_bs, vect_ns, vect_missing_algos)
  df[, "missing"] <- missing_vector
  
  col_perf <- c("M", "sup_norm", "MSE", "bias", "var")
  df["J_bs = 5", col_perf] <- perf_J$avg_perf_5_bs
  df["J_bs = 7", col_perf] <- perf_J$avg_perf_7_bs
  df["J_bs = 11", col_perf] <- perf_J$avg_perf_11_bs
  df["J_bs = 19", col_perf] <- perf_J$avg_perf_19_bs
  df["J_bs = 35", col_perf] <- perf_J$avg_perf_35_bs
  
  df["J_ns = 3", col_perf] <- perf_J$avg_perf_3_ns
  df["J_ns = 5", col_perf] <- perf_J$avg_perf_5_ns
  df["J_ns = 9", col_perf] <- perf_J$avg_perf_9_ns
  df["J_ns = 17", col_perf] <- perf_J$avg_perf_17_ns
  df["J_ns = 33", col_perf] <- perf_J$avg_perf_33_ns
  
  df["CVM_bs", col_perf] <- perf_algos$measures_CVM_bs
  df["CVM_ns", col_perf] <- perf_algos$measures_CVM_ns
  df["CVMSE_bs", col_perf] <- perf_algos$measures_CVMSE_bs
  df["CVMSE_ns", col_perf] <- perf_algos$measures_CVMSE_ns
  df["lepski_bs", col_perf] <- perf_algos$measures_lepski_bs
  df["lepski_ns", col_perf] <- perf_algos$measures_lepski_ns
  df["lepskiboot_bs", col_perf] <- perf_algos$measures_lepski_bs
  df["lepskiboot_ns", col_perf] <- perf_algos$measures_lepski_ns
  
  return(list(df_perf = df, for_curves = res_1_filter, 
              values_bs_allJ = values_bs_allJ,
              values_ns_allJ = values_ns_allJ))
  
  # return(list(res_filtered = res_1_filter, perf_algos = perf_algos,
  #             n_CVM_bs = n_CVM_bs, n_CVM_ns = n_CVM_ns,
  #             n_CVMSE_bs = n_CVMSE_bs, n_CVMSE_ns = n_CVMSE_ns,
  #             n_lepski_bs = n_lepski_bs, n_lepski_ns = n_lepski_ns,
  #             n_lepskiboot_bs = n_lepskiboot_bs, n_lepskiboot_ns = n_lepskiboot_ns,
  #             missing_iter_J = missing_iter,
  #             perf_J = perf_J))
}

#### test on a simulation given ####
n_MC = 500
degree = 3
p_train = 0.5
n_boot = 100
rhozw = 0.9
rhouv = 0.5
case = 2
n_values = 1000

J_bs <- c(5, 7, 11, 19, 35)
J_ns <- c(3, 5, 9, 17, 33)

n_eval = 100

g_on_x = g_sim_3(seq(-2, 2, length.out = 100), 2)


res <- create_df_measures(n_MC, degree, p_train, n_boot, rhozw, rhouv, case, n_values,g_on_x)
df <- res$df_perf
for_curves <- res$for_curves

#save(for_curves, file = "for_curves_test")
df
