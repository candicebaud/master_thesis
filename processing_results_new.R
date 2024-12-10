#### Packages ####
#setwd("C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/github/master_thesis")
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

#### Load results ####
#load("opt_5_degree3_ptrain0.5_nboot100_rhozw0.9_rhouv0.5_case2_n1000.R")
file_list <- list.files(path = "C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/github/master_thesis",
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

#### Create the list containing results to analyze ####
J_bs <- c(5, 7, 11, 19, 35)
J_ns <- c(3, 5, 9, 17, 33)

n_eval = 100

create_list_to_analyze <- function(n_MC, degree, p_train, n_boot, rhozw, rhouv, case, n_values){
  simu_name = paste("opt_", n_MC , "_degree", degree, "_ptrain", p_train, "_nboot",
                    n_boot, "_rhozw" , rhozw,"_rhouv", rhouv , "_case", case,
                    "_n", n_values ,sep = "") 
  
  number_of_simu = length(all_data)
  n_MC = length(all_data[[1]])
  
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
  
  for (n in 1:n_MC){
    all_lists$list_J_opt_CV_M_bs[n] <- all_data[[simu_name]][[n]]$J_opt_CV_M_bs
    all_lists$list_J_opt_CV_M_ns[n] <- all_data[[simu_name]][[n]]$J_opt_CV_M_ns
    all_lists$list_J_opt_CVMSE_bs[n] <- all_data[[simu_name]][[n]]$J_opt_CVMSE_bs
    all_lists$list_J_opt_CVMSE_ns[n] <- all_data[[simu_name]][[n]]$J_opt_CVMSE_ns
    all_lists$list_J_opt_lepski_bs[n] <- all_data[[simu_name]][[n]]$J_opt_lepski_bs
    all_lists$list_J_opt_lepski_ns[n] <- all_data[[simu_name]][[n]]$J_opt_lepski_ns
    all_lists$list_J_opt_lepskiboot_bs[n] <- all_data[[simu_name]][[n]]$J_opt_lepskiboot_bs
    all_lists$list_J_opt_lepskiboot_ns[n] <- all_data[[simu_name]][[n]]$J_opt_lepskiboot_ns
    all_lists$list_gamma_bs[[n]] <- all_data[[simu_name]][[n]]$gamma_bs
    all_lists$list_gamma_ns[[n]] <- all_data[[simu_name]][[n]]$gamma_ns
    all_lists$list_g_hat_J_bs[[n]] <- all_data[[simu_name]][[n]]$g_hat_J_bs
    all_lists$list_g_hat_J_ns[[n]] <- all_data[[simu_name]][[n]]$g_hat_J_ns
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
  names_to_remove <- c("list_gamma_bs", "list_gamma_ns", "list_g_hat_J_bs", "list_g_hat_J_ns")
  all_lists <- all_lists[setdiff(names(all_lists), names_to_remove)]
  
  return(all_lists)
  
}



#### Test ####
n_MC = 5
degree = 3
p_train = 0.5
n_boot = 100
rhozw = 0.9
rhouv = 0.5
case = 2
n_values = 1000
res2 <- create_list_to_analyze(n_MC, degree, p_train, n_boot, rhozw, rhouv, case, n_values)


#load("C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/simulation_results/1/simu_441463_complete_allJ_benchmark/MC2000_fixedJ_J4_degree3_rhozw0.9_rhouv0.8_case2_n400.R")


#### Adapt the previous functions to new data structure (compute perf et compute new_MC dans ce cas là )####ù


#faire les perf et les missing runs pour chaque méthode 