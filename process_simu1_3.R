#### Process results ####
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

setwd("C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/simulation_results/5_final")

#### Import data ####
# file_list <- list.files(path = "C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/simulation_results/4_final",
#                         pattern = "^opt_.*\\.R$", full.names = FALSE)
# 
# all_data <- lapply(file_list, function(filepath) {
#   # Extract the filename without the full path
#   filename <- basename(filepath)
#   name <- str_remove(filename, "\\.R$")  # Remove the ".R" extension
#   
#   # Load the data from the file
#   load(filepath)  # Assumes "results" object is loaded
#   
#   # Return the results, named by the filename
#   list(name = name, data = results)
# })
# 
# all_data <- setNames(lapply(all_data, `[[`, "data"), sapply(all_data, `[[`, "name"))


#### Missing number ####
missing_number <- function(simu_est, n_MC, degree, p_train, n_boot, rhozw, rhouv, case, n_values){
  # simu_name = paste("opt_", n_MC , "_degree", degree, "_ptrain", p_train, "_nboot",
  #                   n_boot, "_rhozw" , rhozw,"_rhouv", rhouv , "_case", case,
  #                   "_n", n_values ,sep = "") 
  
  is_null <- sapply(simu_est, is.null)
  n_null <- sum(is_null) #first the buggs so it did not even run : so did not run for any method
  simu_est <- simu_est[!is_null] #remove them
  n_MC = length(simu_est)
  
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
    list_gamma_bs <- simu_est[[n]]$list_gamma_bs
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
    
    list_gamma_ns <- simu_est[[n]]$list_gamma_ns
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
              n_gamma_ns_33 = n_gamma_ns_33, is_null_list = is_null)) #on renvoie la liste des null pour ensuite ajouter la data générée
}

#### all lists definition ####
reshape_simul_all <- function(simul_all){#simul_all déjà filtré les null donc ok normalement
  n_MC = length(simul_all)
  n_values = length(simul_all[[1]]$Y)
  simul_all_sorted <- list(length=3)
  for (i in 1:2){
    simul_all_sorted$Y <- matrix(0, n_MC, n_values)
    simul_all_sorted$W <- matrix(0, n_MC, n_values)
    simul_all_sorted$Z <- matrix(0, n_MC, n_values)
  }
  for (n in 1:n_MC){
    simul_all_sorted$Y[n, ] <- simul_all[[n]]$Y
    simul_all_sorted$W[n, ] <- simul_all[[n]]$W 
    simul_all_sorted$Z[n, ] <- simul_all[[n]]$Z 
  }
  return(simul_all_sorted)
}


create_list_to_analyze <- function(simu_est, n_MC, degree, p_train, n_boot, rhozw, rhouv, case, n_values, simul_all, is_null){
  # simu_name = paste("opt_", n_MC , "_degree", degree, "_ptrain", p_train, "_nboot",
  #                   n_boot, "_rhozw" , rhozw,"_rhouv", rhouv , "_case", case,
  #                   "_n", n_values ,sep = "") 
  #is_null <- sapply(all_data[[simu_name]], is.null) #length 2000 attention, faudra garder les indices et checker avec la data
  simu_est <- simu_est[!is_null]
  
  #simul_all reshaped
  simul_all <- simul_all[!is_null]
  simul_all_reshaped <- reshape_simul_all(simul_all)
  
  #new length
  n_MC = length(simu_est)
  
  #create list with all results
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
  all_lists$list_matrix_g_hat_J_bs <- vector("list", length = 5) #one matrix per parameter that will form a list
  all_lists$list_matrix_g_hat_J_ns <- vector("list", length = 5)
  all_lists$matrix_W <- simul_all_reshaped$W
  all_lists$matrix_Y <-  simul_all_reshaped$Y
  all_lists$matrix_Z <-  simul_all_reshaped$Z
  
  #à J fixé, on peut faire une liste de matrices par méthode 
  all_lists$list_g_hat_Z_bs_J_fixed <- vector("list", length = 5)
  all_lists$list_g_hat_Z_ns_J_fixed <- vector("list", length = 5)
  
  for (n in 1:n_MC){
    all_lists$list_J_opt_CV_M_bs[n] <- simu_est[[n]]$J_opt_CV_M_bs
    all_lists$list_J_opt_CV_M_ns[n] <- simu_est[[n]]$J_opt_CV_M_ns
    all_lists$list_J_opt_CVMSE_bs[n] <- simu_est[[n]]$J_opt_CVMSE_bs
    all_lists$list_J_opt_CVMSE_ns[n] <- simu_est[[n]]$J_opt_CVMSE_ns
    all_lists$list_J_opt_lepski_bs[n] <- simu_est[[n]]$J_opt_lepski_bs
    all_lists$list_J_opt_lepski_ns[n] <- simu_est[[n]]$J_opt_lepski_ns
    all_lists$list_J_opt_lepskiboot_bs[n] <- simu_est[[n]]$J_opt_lepskiboot_bs
    all_lists$list_J_opt_lepskiboot_ns[n] <- simu_est[[n]]$J_opt_lepskiboot_ns
    all_lists$list_gamma_bs[[n]] <- simu_est[[n]]$list_gamma_bs
    all_lists$list_gamma_ns[[n]] <- simu_est[[n]]$list_gamma_ns
    all_lists$list_g_hat_J_bs[[n]] <- simu_est[[n]]$list_g_hat_J_bs
    all_lists$list_g_hat_J_ns[[n]] <- simu_est[[n]]$list_g_hat_J_ns
  }
  
  #save the generated data
  for (i in 1:5){
    all_lists$list_matrix_g_hat_J_bs[[i]] <- matrix(0, nrow = n_MC, ncol = n_eval)
    all_lists$list_matrix_g_hat_J_ns[[i]] <- matrix(0, nrow = n_MC, ncol = n_eval)
    
    all_lists$list_g_hat_Z_bs_J_fixed[[i]] <- matrix(0, nrow = n_MC, ncol = 1000)
    all_lists$list_g_hat_Z_ns_J_fixed[[i]] <- matrix(0, nrow = n_MC, ncol = 1000)
    
    for (n in 1:n_MC){
      all_lists$list_matrix_g_hat_J_bs[[i]][n,] <- simu_est[[n]]$list_g_hat_J_bs[[i]][,1]
      all_lists$list_matrix_g_hat_J_ns[[i]][n,] <- simu_est[[n]]$list_g_hat_J_ns[[i]][,1]
      
      # #créer base de splines 
      # Z_n <- all_lists$matrix_Z[n,]
      # P_bs <- create_dyadic_P_splines_bs(Z_n, Z_n, as.numeric(J_bs[i]), 3)
      # P_ns <- create_dyadic_P_splines_ns(Z_n, Z_n, as.numeric(J_ns[i]), 3)
      # 
      # #récupérer le gamma et faire le produit matriciel
      # gamma_bs_J <- all_lists$list_gamma_bs[[n]][[i]]
      # gamma_ns_J <- all_lists$list_gamma_ns[[n]][[i]]
      # 
      # #stocker la pred
      # all_lists$list_g_hat_Z_bs_J_fixed[[i]][n,] <- P_bs%*%gamma_bs_J
      # all_lists$list_g_hat_Z_ns_J_fixed[[i]][n,] <- P_ns%*%gamma_ns_J
      }
  }
  
  #create the lists for each method of values on x_eval obtained
  all_lists$list_est_values_CV_M_bs <- matrix(0, nrow = n_MC, ncol = n_eval)
  all_lists$list_est_values_CV_M_ns <- matrix(0, nrow = n_MC, ncol = n_eval)
  all_lists$list_est_values_CVMSE_bs <- matrix(0, nrow = n_MC, ncol = n_eval)
  all_lists$list_est_values_CVMSE_ns <- matrix(0, nrow = n_MC, ncol = n_eval)
  all_lists$list_est_values_lepski_bs <- matrix(0, nrow = n_MC, ncol = n_eval)
  all_lists$list_est_values_lepski_ns <- matrix(0, nrow = n_MC, ncol = n_eval)
  all_lists$list_est_values_lepskiboot_bs <- matrix(0, nrow = n_MC, ncol = n_eval)
  all_lists$list_est_values_lepskiboot_ns <- matrix(0, nrow = n_MC, ncol = n_eval)
  
  all_lists$list_gamma_CV_M_bs <- vector("list", length = n_MC)
  all_lists$list_gamma_CVMSE_bs <- vector("list", length = n_MC)
  all_lists$list_gamma_lepski_bs <- vector("list", length = n_MC)
  all_lists$list_gamma_lepskiboot_bs <- vector("list", length = n_MC)
  all_lists$list_gamma_CV_M_ns <- vector("list", length = n_MC)
  all_lists$list_gamma_CVMSE_ns <- vector("list", length = n_MC)
  all_lists$list_gamma_lepski_ns <- vector("list", length = n_MC)
  all_lists$list_gamma_lepskiboot_ns <- vector("list", length = n_MC)
  
  
  # fill the matrix with the estimated values
  for (n in 1:n_MC){
    #bs
    index = which(J_bs == all_lists$list_J_opt_CV_M_bs[n])
    if (length(index)>0){
      all_lists$list_est_values_CV_M_bs[n,] <- all_lists$list_g_hat_J_bs[[n]][index][[1]]
      all_lists$list_gamma_CV_M_bs[[n]] <- all_lists$list_gamma_bs[[n]][index][[1]]}
      
    index = which(J_bs == all_lists$list_J_opt_CVMSE_bs[n])
    if (length(index)>0){
      all_lists$list_est_values_CVMSE_bs[n,] <- all_lists$list_g_hat_J_bs[[n]][index][[1]]
      all_lists$list_gamma_CVMSE_bs[[n]] <- all_lists$list_gamma_bs[[n]][index][[1]] }
    
    index = which(J_bs == all_lists$list_J_opt_lepski_bs[n])
    if (length(index)>0){
      all_lists$list_est_values_lepski_bs[n,] <- all_lists$list_g_hat_J_bs[[n]][index][[1]]
      all_lists$list_gamma_lepski_bs[[n]] <- all_lists$list_gamma_bs[[n]][index][[1]]}
    
    index = which(J_bs == all_lists$list_J_opt_lepskiboot_bs[n])
    if (length(index)>0){
      all_lists$list_est_values_lepskiboot_bs[n,] <- all_lists$list_g_hat_J_bs[[n]][index][[1]]
      all_lists$list_gamma_lepskiboot_bs[[n]] <- all_lists$list_gamma_bs[[n]][index][[1]]}
    
    index = which(J_ns == all_lists$list_J_opt_CV_M_ns[n])
    if (length(index)>0){
      all_lists$list_est_values_CV_M_ns[n,] <- all_lists$list_g_hat_J_ns[[n]][index][[1]]
      all_lists$list_gamma_CV_M_ns[[n]] <- all_lists$list_gamma_ns[[n]][index][[1]]}
      
    index = which(J_ns == all_lists$list_J_opt_CVMSE_ns[n])
    if (length(index)>0){
      all_lists$list_est_values_CVMSE_ns[n,] <- all_lists$list_g_hat_J_ns[[n]][index][[1]]
      all_lists$list_gamma_CVMSE_ns[[n]] <- all_lists$list_gamma_ns[[n]][index][[1]]}
      
    index = which(J_ns == all_lists$list_J_opt_lepski_ns[n])
    if (length(index)>0){
      all_lists$list_est_values_lepski_ns[n,] <- all_lists$list_g_hat_J_ns[[n]][index][[1]]
      all_lists$list_gamma_lepski_ns[[n]] <- all_lists$list_gamma_ns[[n]][index][[1]]}  
      
    index = which(J_ns == all_lists$list_J_opt_lepskiboot_ns[n])
    if (length(index)>0){
      all_lists$list_est_values_lepskiboot_ns[n,] <- all_lists$list_g_hat_J_ns[[n]][index][[1]]
      all_lists$list_gamma_lepskiboot_ns[[n]] <- all_lists$list_gamma_ns[[n]][index][[1]]}
  }
  
  return(all_lists)
}


#### Create lists for each method of estimation
filter_res <- function(res_to_analyze){
  #CVM bs
  res_CVM_bs <- list()
  zero_indices <- which(res_to_analyze$list_J_opt_CV_M_bs == 0) #les J qui valent 0 
  n_CVM_bs <- length(zero_indices)
  if (n_CVM_bs>0){
    res_CVM_bs$list_J_opt_CV_M_bs <- res_to_analyze$list_J_opt_CV_M_bs[-zero_indices]
    res_CVM_bs$list_est_values_CV_M_bs <- res_to_analyze$list_est_values_CV_M_bs[-zero_indices,]
    res_CVM_bs$matrix_W <- res_to_analyze$matrix_W[-zero_indices,]
    res_CVM_bs$matrix_Z <- res_to_analyze$matrix_Z[-zero_indices,]
    res_CVM_bs$matrix_Y <- res_to_analyze$matrix_Y[-zero_indices,]
    res_CVM_bs$list_gamma <- res_to_analyze$list_gamma_CV_M_bs[-zero_indices]
  }
  else{
    res_CVM_bs$list_J_opt_CV_M_bs <- res_to_analyze$list_J_opt_CV_M_bs
    res_CVM_bs$list_est_values_CV_M_bs <- res_to_analyze$list_est_values_CV_M_bs
    res_CVM_bs$matrix_W <- res_to_analyze$matrix_W
    res_CVM_bs$matrix_Z <- res_to_analyze$matrix_Z
    res_CVM_bs$matrix_Y <- res_to_analyze$matrix_Y
    res_CVM_bs$list_gamma <- res_to_analyze$list_gamma_CV_M_bs
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
    res_CVM_ns$list_gamma <- res_to_analyze$list_gamma_CV_M_ns[-zero_indices]
  }
  else{
    res_CVM_ns$list_J_opt_CV_M_ns <- res_to_analyze$list_J_opt_CV_M_ns
    res_CVM_ns$list_est_values_CV_M_ns <- res_to_analyze$list_est_values_CV_M_ns
    res_CVM_ns$matrix_W <- res_to_analyze$matrix_W
    res_CVM_ns$matrix_Z <- res_to_analyze$matrix_Z
    res_CVM_ns$matrix_Y <- res_to_analyze$matrix_Y
    res_CVM_ns$list_gamma <- res_to_analyze$list_gamma_CV_M_ns
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
    res_CVMSE_bs$list_gamma <- res_to_analyze$list_gamma_CVMSE_bs[-zero_indices]
  }
  else{
    res_CVMSE_bs$list_J_opt_CVMSE_bs <- res_to_analyze$list_J_opt_CVMSE_bs
    res_CVMSE_bs$list_est_values_CVMSE_bs <- res_to_analyze$list_est_values_CVMSE_bs
    res_CVMSE_bs$matrix_W <- res_to_analyze$matrix_W
    res_CVMSE_bs$matrix_Z <- res_to_analyze$matrix_Z
    res_CVMSE_bs$matrix_Y <- res_to_analyze$matrix_Y
    res_CVMSE_bs$list_gamma <- res_to_analyze$list_gamma_CVMSE_bs
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
    res_CVMSE_ns$list_gamma <- res_to_analyze$list_gamma_CVMSE_ns[-zero_indices]
  }
  else{
    res_CVMSE_ns$list_J_opt_CVMSE_ns <- res_to_analyze$list_J_opt_CVMSE_ns
    res_CVMSE_ns$list_est_values_CVMSE_ns <- res_to_analyze$list_est_values_CVMSE_ns
    res_CVMSE_ns$matrix_W <- res_to_analyze$matrix_W
    res_CVMSE_ns$matrix_Z <- res_to_analyze$matrix_Z
    res_CVMSE_ns$matrix_Y <- res_to_analyze$matrix_Y
    res_CVMSE_ns$list_gamma <- res_to_analyze$list_gamma_CVMSE_ns
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
    res_lepski_bs$list_gamma <- res_to_analyze$list_gamma_lepski_bs[-zero_indices]
  }
  else{
    res_lepski_bs$list_J_opt_lepski_bs <- res_to_analyze$list_J_opt_lepski_bs
    res_lepski_bs$list_est_values_lepski_bs <- res_to_analyze$list_est_values_lepski_bs
    res_lepski_bs$matrix_W <- res_to_analyze$matrix_W
    res_lepski_bs$matrix_Z <- res_to_analyze$matrix_Z
    res_lepski_bs$matrix_Y <- res_to_analyze$matrix_Y
    res_lepski_bs$list_gamma <- res_to_analyze$list_gamma_lepski_bs
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
    res_lepski_ns$list_gamma <- res_to_analyze$list_gamma_lepski_ns[-zero_indices]
  }
  else{
    res_lepski_ns$list_J_opt_lepski_ns <- res_to_analyze$list_J_opt_lepski_ns
    res_lepski_ns$list_est_values_lepski_ns <- res_to_analyze$list_est_values_lepski_ns
    res_lepski_ns$matrix_W <- res_to_analyze$matrix_W
    res_lepski_ns$matrix_Z <- res_to_analyze$matrix_Z
    res_lepski_ns$matrix_Y <- res_to_analyze$matrix_Y
    res_lepski_ns$list_gamma <- res_to_analyze$list_gamma_lepski_ns
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
    res_lepskiboot_bs$list_gamma <- res_to_analyze$list_gamma_lepskiboot_bs[-zero_indices]
  }
  else{
    res_lepskiboot_bs$list_J_opt_lepskiboot_bs <- res_to_analyze$list_J_opt_lepskiboot_bs
    res_lepskiboot_bs$list_est_values_lepskiboot_bs <- res_to_analyze$list_est_values_lepskiboot_bs
    res_lepskiboot_bs$matrix_W <- res_to_analyze$matrix_W
    res_lepskiboot_bs$matrix_Z <- res_to_analyze$matrix_Z
    res_lepskiboot_bs$matrix_Y <- res_to_analyze$matrix_Y
    res_lepskiboot_bs$list_gamma <- res_to_analyze$list_gamma_lepskiboot_bs
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
    res_lepskiboot_ns$list_gamma <- res_to_analyze$list_gamma_lepskiboot_ns[-zero_indices]
  }
  else{
    res_lepskiboot_ns$list_J_opt_lepskiboot_ns <- res_to_analyze$list_J_opt_lepskiboot_ns
    res_lepskiboot_ns$list_est_values_lepskiboot_ns <- res_to_analyze$list_est_values_lepskiboot_ns
    res_lepskiboot_ns$matrix_W <- res_to_analyze$matrix_W
    res_lepskiboot_ns$matrix_Z <- res_to_analyze$matrix_Z
    res_lepskiboot_ns$matrix_Y <- res_to_analyze$matrix_Y
    res_lepskiboot_ns$list_gamma <- res_to_analyze$list_gamma_lepskiboot_ns
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


#### Compute perf for estimation methods ####
compute_all_perf <- function(mat_g_hat_on_x, g_0_on_x, data_algo){
  n_eval = length(g_0_on_x)
  n_MC = nrow(mat_g_hat_on_x)
  
  if (n_MC>0){
    avg <- colMeans(mat_g_hat_on_x)
    
    MSE = sum(rowSums((mat_g_hat_on_x - matrix(g_0_on_x, nrow = n_MC, ncol = n_eval, byrow = TRUE))^2)) / (n_MC * n_eval)
    
    var <- sum(rowSums((mat_g_hat_on_x - matrix(avg, nrow = n_MC, ncol = n_eval, byrow = TRUE))^2)) / (n_MC * n_eval)
    
    bias <- sum((g_0_on_x - avg)^2) / n_eval
    
    sup_norm_vect <- apply(abs(mat_g_hat_on_x - matrix(g_0_on_x, nrow = n_MC, ncol = n_eval, byrow = TRUE)), 1, max)
    sup_norm <- mean(sup_norm_vect)
    
    #M 
    M_vect <- rep(0, n_MC)
    matrix_Y <- data_algo$matrix_Y
    matrix_Z <- data_algo$matrix_Z
    matrix_W <- data_algo$matrix_W
    
    #pour calculer M, il faut que je calcule g_hat_Z soit pour les J soit pour la m"thode : à rajouter dans chaque algo pour avoir une liste 
    #et ensuite j'utilise la formule du papier de lapenta
    
    # Omega_list <- lapply(1:n_MC, function(n) create_W(matrix_W[n, ]))
    # 
    # residuals <- matrix_Y - mat_g_hat_on_x
    # M_vect <- sapply(1:n_MC, function(n) {
    #   residual <- residuals[n, ]
    #   Omega <- Omega_list[[n]]
    #   t(residual) %*% Omega %*% residual / (n_eval^2)
    # })
    M = mean(M_vect)
    
    #list(MSE = MSE, var = var, bias = bias, supnorm = sup_norm, M = M)
    return(c(M, sup_norm, MSE, bias, var))}
  else{
    return(rep(999, 5))
  }
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


##### Compute perf for J fixed ####
filter_matrix_pred <- function(matrix_pred) {
  # Check which rows have a sum not equal to 0
  non_zero_rows <- apply(matrix_pred, 1, function(row) sum(row) != 0)
  
  # Subset the matrix to include only the non-zero rows
  new_matrix <- matrix_pred[non_zero_rows, , drop = FALSE]
  
  return(list(matrix = new_matrix, non_zero_rows = non_zero_rows))
}

compute_perf_J <- function(g_on_x, all_lists){ #TO DO
  value_bs_all_J <- all_lists$list_matrix_g_hat_J_bs
  value_ns_all_J <- all_lists$list_matrix_g_hat_J_ns
  
  #initialize the lists that will contain results specific to each J 
  #dans chaque j'aurai un élément : g_hat_on_x en matrice, Y et W en matrice
  list_5_bs <- list(length = 3)
  list_7_bs <- list(length = 3)
  list_11_bs <- list(length = 3)
  list_19_bs <- list(length = 3)
  list_35_bs <- list(length = 3)
  
  list_3_ns <- list(length = 3)
  list_5_ns <- list(length = 3)
  list_9_ns <- list(length = 3)
  list_17_ns <- list(length = 3)
  list_33_ns <- list(length = 3)
  
  #filter the indices for bad simulations and fill the sublists
  for (i in 1:5){
    filter_bs_all_J <- filter_matrix_pred(value_bs_all_J[[i]])
    filter_ns_all_J <- filter_matrix_pred(value_ns_all_J[[i]])
    if (i == 1){
      list_5_bs$mat_g_hat_on_x <- filter_bs_all_J$matrix
      non_zero_rows <- filter_bs_all_J$non_zero_rows
      list_5_bs$matrix_Y <- all_lists$matrix_Y[non_zero_rows, , drop = FALSE]
      list_5_bs$matrix_W <- all_lists$matrix_W[non_zero_rows, , drop = FALSE]
      
      list_3_ns$mat_g_hat_on_x <- filter_ns_all_J$matrix
      non_zero_rows <- filter_ns_all_J$non_zero_rows
      list_3_ns$matrix_Y <- all_lists$matrix_Y[non_zero_rows, , drop = FALSE]
      list_3_ns$matrix_W <- all_lists$matrix_W[non_zero_rows, , drop = FALSE]
    }
    
    if (i == 2){
      list_7_bs$mat_g_hat_on_x <- filter_bs_all_J$matrix
      non_zero_rows <- filter_bs_all_J$non_zero_rows
      list_7_bs$matrix_Y <- all_lists$matrix_Y[non_zero_rows, , drop = FALSE]
      list_7_bs$matrix_W <- all_lists$matrix_W[non_zero_rows, , drop = FALSE]
      
      list_5_ns$mat_g_hat_on_x <- filter_ns_all_J$matrix
      non_zero_rows <- filter_ns_all_J$non_zero_rows
      list_5_ns$matrix_Y <- all_lists$matrix_Y[non_zero_rows, , drop = FALSE]
      list_5_ns$matrix_W <- all_lists$matrix_W[non_zero_rows, , drop = FALSE]
    }
    
    if (i == 3){
      list_11_bs$mat_g_hat_on_x <- filter_bs_all_J$matrix
      non_zero_rows <- filter_bs_all_J$non_zero_rows
      list_11_bs$matrix_Y <- all_lists$matrix_Y[non_zero_rows, , drop = FALSE]
      list_11_bs$matrix_W <- all_lists$matrix_W[non_zero_rows, , drop = FALSE]
      
      list_9_ns$mat_g_hat_on_x <- filter_ns_all_J$matrix
      non_zero_rows <- filter_ns_all_J$non_zero_rows
      list_9_ns$matrix_Y <- all_lists$matrix_Y[non_zero_rows, , drop = FALSE]
      list_9_ns$matrix_W <- all_lists$matrix_W[non_zero_rows, , drop = FALSE]
    }
    
    if (i == 4){
      list_19_bs$mat_g_hat_on_x <- filter_bs_all_J$matrix
      non_zero_rows <- filter_bs_all_J$non_zero_rows
      list_19_bs$matrix_Y <- all_lists$matrix_Y[non_zero_rows, , drop = FALSE]
      list_19_bs$matrix_W <- all_lists$matrix_W[non_zero_rows, , drop = FALSE]
      
      list_17_ns$mat_g_hat_on_x <- filter_ns_all_J$matrix
      non_zero_rows <- filter_ns_all_J$non_zero_rows
      list_17_ns$matrix_Y <- all_lists$matrix_Y[non_zero_rows, , drop = FALSE]
      list_17_ns$matrix_W <- all_lists$matrix_W[non_zero_rows, , drop = FALSE]
    }
    
    if (i == 5){
      list_35_bs$mat_g_hat_on_x <- filter_bs_all_J$matrix
      non_zero_rows <- filter_bs_all_J$non_zero_rows
      list_35_bs$matrix_Y <- all_lists$matrix_Y[non_zero_rows, , drop = FALSE]
      list_35_bs$matrix_W <- all_lists$matrix_W[non_zero_rows, , drop = FALSE]
      
      list_33_ns$mat_g_hat_on_x <- filter_ns_all_J$matrix
      non_zero_rows <- filter_ns_all_J$non_zero_rows
      list_33_ns$matrix_Y <- all_lists$matrix_Y[non_zero_rows, , drop = FALSE]
      list_33_ns$matrix_W <- all_lists$matrix_W[non_zero_rows, , drop = FALSE]
    }
  }
  
  to_return <- list()
  to_return$avg_perf_5_bs <- compute_all_perf(list_5_bs$mat_g_hat_on_x, 
                                              g_on_x, list_5_bs)
  to_return$avg_perf_7_bs <- compute_all_perf(list_7_bs$mat_g_hat_on_x, 
                                              g_on_x, list_7_bs)
  to_return$avg_perf_11_bs <- compute_all_perf(list_11_bs$mat_g_hat_on_x, 
                                               g_on_x, list_11_bs)
  to_return$avg_perf_19_bs <- compute_all_perf(list_19_bs$mat_g_hat_on_x, 
                                               g_on_x, list_19_bs)
  to_return$avg_perf_35_bs <- compute_all_perf(list_35_bs$mat_g_hat_on_x, 
                                               g_on_x, list_35_bs)
  
  to_return$avg_perf_3_ns <- compute_all_perf(list_3_ns$mat_g_hat_on_x, 
                                              g_on_x, list_3_ns)
  to_return$avg_perf_5_ns <- compute_all_perf(list_5_ns$mat_g_hat_on_x, 
                                              g_on_x, list_5_ns)
  to_return$avg_perf_9_ns <- compute_all_perf(list_9_ns$mat_g_hat_on_x, 
                                              g_on_x, list_9_ns)
  to_return$avg_perf_17_ns <- compute_all_perf(list_17_ns$mat_g_hat_on_x, 
                                               g_on_x, list_17_ns)
  to_return$avg_perf_33_ns <- compute_all_perf(list_33_ns$mat_g_hat_on_x, 
                                               g_on_x, list_33_ns)
  

  return(to_return)
  
}

#### Data process all function ####
create_df_measures <- function(simu_est, n_MC, degree, p_train, n_boot, rhozw, rhouv, case, n_values, g_on_x){
  #load simulation data 
  setwd("C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/simulation_results/5_final")
  data_name = paste("data_", n_MC, "_rhozw" , rhozw,"_rhouv", rhouv , "_case", case, "_n", n_values, ".R" ,sep = "")
  simul_all <- get(load(data_name))
  
  #missing data 
  missing_iter <- missing_number(simu_est, n_MC, degree, p_train, n_boot, rhozw, rhouv, case, n_values)
  n_null = missing_iter$n_null #not even computed
  is_null = missing_iter$is_null
  
  #algorithms, filter 
  all_lists <- create_list_to_analyze(simu_est, n_MC, degree, p_train, n_boot, rhozw, rhouv, case, n_values, simul_all, is_null)
  lists_method <- filter_res(all_lists)
  
  
  #perf algos
  perf_algos <- compute_perf(lists_method, g_on_x)
  
  
  #perf J fixed
  perf_J <- compute_perf_J(g_on_x, all_lists)
  
  
  # missing iterations on the algos
  n_CVM_bs = lists_method$n_CVM_bs + n_null
  n_CVM_ns = lists_method$n_CVM_ns + n_null
  n_CVMSE_bs = lists_method$n_CVMSE_bs + n_null
  n_CVMSE_ns = lists_method$n_CVMSE_ns + n_null
  n_lepski_bs = lists_method$n_lepski_bs + n_null
  n_lepski_ns = lists_method$n_lepski_ns + n_null
  n_lepskiboot_bs = lists_method$n_lepskiboot_bs + n_null
  n_lepskiboot_ns = lists_method$n_lepskiboot_ns + n_null
  
  
  vect_missing_algos <- c(n_CVM_bs, n_CVM_ns, n_CVMSE_bs, n_CVMSE_ns,
                          n_lepski_bs, n_lepski_ns, n_lepskiboot_bs, n_lepskiboot_ns)
  
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
  df["lepskiboot_bs", col_perf] <- perf_algos$measures_lepskiboot_bs
  df["lepskiboot_ns", col_perf] <- perf_algos$measures_lepskiboot_ns
  
  return(list(df_perf = df, for_curves = lists_method,
              all_lists = all_lists
              #values_bs_allJ = values_bs_allJ,
              #values_ns_allJ = values_ns_allJ
  ))
}


#### Analyze results now ####
#simulation 1
n_MC = 2000
degree = 3
p_train = 0.5
n_boot = 100
rhozw = 0.9
rhouv = 0.5
case = 2
n_values = 1000

n_eval = 100

J_bs <- c(4, 5, 7, 11, 19)
J_ns <- c(2, 3, 5, 9, 17)

g_on_x = g_sim_3(seq(-2, 2, length.out = 100), 2)

simu_name = paste("opt_", n_MC , "_degree", degree, "_ptrain", p_train, "_nboot",
                  n_boot, "_rhozw" , rhozw,"_rhouv", rhouv , "_case", case,
                  "_n", n_values, ".R", sep = "") 
simu_est <- get(load(simu_name))

res_1 <- create_df_measures(simu_est, n_MC, degree, p_train, n_boot, rhozw, rhouv, case, n_values,g_on_x)
#save(res_1, file = "res_1")
#library(xtable)
xtable(res_1$df_perf, caption = "Simulation 1" )


#simulation 3
n_MC = 2000
degree = 3
p_train = 0.5
n_boot = 100
rhozw = 0.9
rhouv = 0.5
case = 3
n_values = 1000

n_eval = 100

J_bs <- c(4, 5, 7, 11, 19)
J_ns <- c(2, 3, 5, 9, 17)

g_on_x = g_sim_3(seq(-2, 2, length.out = 100), 3)

simu_name = paste("opt_", n_MC , "_degree", degree, "_ptrain", p_train, "_nboot",
                  n_boot, "_rhozw" , rhozw,"_rhouv", rhouv , "_case", case,
                  "_n", n_values, ".R", sep = "") 
simu_est <- get(load(simu_name))

res_3 <- create_df_measures(simu_est, n_MC, degree, p_train, n_boot, rhozw, rhouv, case, n_values,g_on_x)
save(res_3, file = "res_3")
#library(xtable)
xtable(res_3$df_perf, caption = "Simulation 3" )


#### Compute M for J fixed ####
compute_things_for_M_J_fixed <- function(all_lists){
  n_MC = length(all_lists$list_J_opt_CV_M_bs)
  new_list <- list()
  new_list$list_g_hat_Z_bs_J_fixed <- all_lists$list_g_hat_Z_bs_J_fixed
  new_list$list_g_hat_Z_ns_J_fixed <- all_lists$list_g_hat_Z_ns_J_fixed
  
  new_list$list_gamma_bs <- all_lists$list_gamma_bs
  new_list$list_gamma_ns <- all_lists$list_gamma_ns
  
  new_list$matrix_Z <- all_lists$matrix_Z
  new_list$matrix_Y <- all_lists$matrix_Y
  new_list$matrix_W <- all_lists$matrix_W
  
  rm(all_lists)
  
  for (i in 1:5){
    print(i)
    for (n in 1:n_MC){
      print(n)
      #créer base de splines
      Z_n <- new_list$matrix_Z[n,]
      P_bs <- create_dyadic_P_splines_bs(Z_n, Z_n, as.numeric(J_bs[i]), 3)
      P_ns <- create_dyadic_P_splines_ns(Z_n, Z_n, as.numeric(J_ns[i]), 3)
      
      #récupérer le gamma et faire le produit matriciel
      gamma_bs_J <- new_list$list_gamma_bs[[n]][[i]]
      gamma_ns_J <- new_list$list_gamma_ns[[n]][[i]]
      
      #stocker la pred
      new_list$list_g_hat_Z_bs_J_fixed[[i]][n,] <- P_bs%*%gamma_bs_J
      new_list$list_g_hat_Z_ns_J_fixed[[i]][n,] <- P_ns%*%gamma_ns_J
    }
  }
  
  return(new_list)
} 

#compute and save to avoid recomputing 
for_M_res_1_J_fixed <- compute_things_for_M_J_fixed(res_1$all_lists)
save(for_M_res_1_J_fixed, file = "for_M_res_1_J_fixed")

for_M_res_3_J_fixed <- compute_things_for_M_J_fixed(res_3$all_lists)
save(for_M_res_3_J_fixed, file = "for_M_res_3_J_fixed")


# n_Omega = nrow(for_M_res_1_J_fixed$matrix_W)
# Omega_list_res_1 <- lapply(1:n_Omega, function(n) create_W(for_M_res_1_J_fixed$matrix_W[n, ]))
# save(Omega_list_res_1, file = "Omega_list_res_1")
# 
# n_Omega = nrow(for_M_res_3_J_fixed$matrix_W)
# Omega_list_res_3 <- lapply(1:n_Omega, function(n) create_W(for_M_res_3_J_fixed$matrix_W[n, ]))
# save(Omega_list_res_3, file = "Omega_list_res_3")

# compute M finally
compute_M <- function(for_M){
  n_MC = length(for_M$list_gamma_bs)
  M_bs <- rep(0, 5)
  M_ns <- rep(0, 5)
  
  #compute M omitting when gamma = 0 
  for (i in 1:5){
    print(i)
    residuals <- for_M$matrix_Y - for_M$list_g_hat_Z_bs_J_fixed[[i]]
    M_vect <- sapply(1:n_MC, function(n) {
      if (sum(for_M$list_gamma_bs[[n]][[i]])>0){
        residual <- residuals[n, ]
        Omega <- create_W(for_M$matrix_W[n,])
        t(residual) %*% Omega %*% residual / (1000^2)
      }
      else{
        NA
      }})
    M_bs[i] = mean(M_vect, na.rm = TRUE)
    
    residuals <- for_M$matrix_Y - for_M$list_g_hat_Z_ns_J_fixed[[i]]
    M_vect <- sapply(1:n_MC, function(n) {
      if (sum(for_M$list_gamma_ns[[n]][[i]])>0){
        residual <- residuals[n, ]
        Omega <- create_W(for_M$matrix_W[n,])
        t(residual) %*% Omega %*% residual / (1000^2)
      }
      else{
        NA
      }})
    M_ns[i] = mean(M_vect, na.rm = TRUE)
  }
  
  return(list(M_bs = M_bs, M_ns = M_ns))
} #returns the values of M for all J in bs and ns cases 


#load("for_M_res_1_J_fixed")
M_val_res_1 <- compute_M(for_M_res_1_J_fixed)
save(M_val_res_1, file = "M_val_res_1")

load("for_M_res_3_J_fixed")
M_val_res_3 <- compute_M(for_M_res_3_J_fixed)
save(M_val_res_3, file = "M_val_res_3")


##### Compute M for methods ####
# on renvoie "for_curves" dans le res, qui correspond à list_methods -> on peut extraire les données pour chaque méthode

compute_M_sub <- function(res_method, bs_bool, list_J){
  #est_values <- res_method[[2]]
  matrix_W <- res_method$matrix_W
  matrix_Y <- res_method$matrix_Y
  matrix_Z <- res_method$matrix_Z
  
  #gamma
  list_gamma <- res_method$list_gamma
  
  n_MC = nrow(res_method$matrix_W)
  
  M_vect <- sapply(1:n_MC, function(n) {
    print(n)
    gamma <- list_gamma[[n]]
    if (sum(gamma)>0){
      Omega <- create_W(matrix_W[n,])
      if (bs_bool ==1){
        P <- create_dyadic_P_splines_bs(matrix_Z[n,], matrix_Z[n,], list_J[n],3)
      }
      else{
        P <- create_dyadic_P_splines_ns(matrix_Z[n,], matrix_Z[n,], list_J[n],3)
      }
      est_values <- P %*% gamma
      residual <- matrix_Y[n,] - est_values
      #residual <- residuals[n, ]
      t(residual) %*% Omega %*% residual / (1000^2)
    }
    else{
      NA
    }})
  return(mean(M_vect, na.rm = TRUE))
}

  
  
compute_M_methods <- function(for_curves){
  M_CVM_bs <- compute_M_sub(for_curves$res_CVM_bs, 1, for_curves$res_CVM_bs$list_J_opt_CV_M_bs)
  print(M_CVM_bs)
  M_CVMSE_bs <- compute_M_sub(for_curves$res_CVMSE_bs, 1, for_curves$res_CVMSE_bs$list_J_opt_CVMSE_bs)
  print(M_CVMSE_bs)
  M_lepski_bs <- compute_M_sub(for_curves$res_lepski_bs, 1, for_curves$res_lepski_bs$list_J_opt_lepski_bs)
  print(M_lepski_bs)
  M_lepskiboot_bs <- compute_M_sub(for_curves$res_lepskiboot_bs, 1, for_curves$res_lepskiboot_bs$list_J_opt_lepskiboot_bs)
  print(M_lepskiboot_bs)
  
  M_CVM_ns <- compute_M_sub(for_curves$res_CVM_ns, 0, for_curves$res_CVM_ns$list_J_opt_CV_M_ns)
  print(M_CVM_ns)
  M_CVMSE_ns <- compute_M_sub(for_curves$res_CVMSE_ns, 0, for_curves$res_CVMSE_ns$list_J_opt_CVMSE_ns)
  print(M_CVMSE_ns)
  M_lepski_ns <- compute_M_sub(for_curves$res_lepski_ns, 0, for_curves$res_lepski_ns$list_J_opt_lepski_ns)
  print(M_lepski_ns)
  M_lepskiboot_ns <- compute_M_sub(for_curves$res_lepskiboot_ns, 0, for_curves$res_lepskiboot_ns$list_J_opt_lepskiboot_ns)
  print(M_lepskiboot_ns)
   
  return(list(M_CVM_bs = M_CVM_bs, M_CVMSE_bs = M_CVMSE_bs,
              M_lepski_bs = M_lepski_bs, M_lepskiboot_bs = M_lepskiboot_bs,
              M_CVM_ns = M_CVM_ns, M_CVMSE_ns = M_CVMSE_ns, 
              M_lepski_ns = M_lepski_ns, M_lepskiboot_ns = M_lepskiboot_ns))
  
  
}

M_methods_res_1 <- compute_M_methods(res_1$for_curves)
save(M_methods_res_1, file = "M_methods_res_1")

M_methods_res_3 <- compute_M_methods(res_3$for_curves)
save(M_methods_res_3, file = "M_methods_res_3")
